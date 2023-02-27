/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const ssao_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common

uniform sampler2D tDepth;
uniform sampler2D tDepthHalf;
uniform sampler2D tDepthQuarter;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform vec3 uSamples[dNSamples];

uniform mat4 uProjection;
uniform mat4 uInvProjection;

uniform float uRadius[dLevels];
uniform float uBias[dLevels];
uniform float uDistanceFactor;
uniform float uMinDistanceFactor;
uniform bool uSolidBackground;

float smootherstep(float edge0, float edge1, float x) {
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

float noise(const in vec2 coords) {
    float a = 12.9898;
    float b = 78.233;
    float c = 43758.5453;
    float dt = dot(coords, vec2(a,b));
    float sn = mod(dt, PI);
    return abs(fract(sin(sn) * c)); // is abs necessary?
}

vec2 getNoiseVec2(const in vec2 coords) {
    return vec2(noise(coords), noise(coords + vec2(PI, 2.71828)));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

float getDepth(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    #ifdef depthTextureSupport
        return texture2D(tDepth, c).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, c));
    #endif
}

float getMappedDepth(const in vec2 coords, const in vec2 selfCoords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    float d = distance(coords, selfCoords);
    #ifdef depthTextureSupport
        if (d > 0.1) {
            return texture2D(tDepthQuarter, c).r;
        } else if (d > 0.05) {
            return texture2D(tDepthHalf, c).r;
        } else {
            return texture2D(tDepth, c).r;
        }
        return texture2D(tDepth, c).r;
    #else
        if (d > 0.0.1) {
            return unpackRGBAToDepth(texture2D(tDepthQuarter, c));
        } else if (d > 0.05) {
            return unpackRGBAToDepth(texture2D(tDepthHalf, c));
        } else {
            return unpackRGBAToDepth(texture2D(tDepth, c));
        }
        return unpackRGBAToDepth(texture2D(tDepth, c));
    #endif
}

vec3 normalFromDepth(const in float depth, const in float depth1, const in float depth2, vec2 offset1, vec2 offset2) {
    vec3 p1 = vec3(offset1, depth1 - depth);
    vec3 p2 = vec3(offset2, depth2 - depth);

    vec3 normal = cross(p1, p2);
    normal.z = -normal.z;

    return normalize(normal);
}

float getPixelSize(const in vec2 coords, const in float depth) {
    vec3 viewPos0 = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
    vec3 viewPos1 = screenSpaceToViewSpace(vec3(coords + vec2(1.0, 0.0) / uTexSize, depth), uInvProjection);
    return distance(viewPos0, viewPos1);
}

// StarCraft II Ambient Occlusion by [Filion and McNaughton 2008]
void main(void) {
    vec2 invTexSize = 1.0 / uTexSize;
    vec2 selfCoords = gl_FragCoord.xy * invTexSize;

    float selfDepth = getDepth(selfCoords);
    vec2 selfPackedDepth = packUnitIntervalToRG(selfDepth);

    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(0.0), selfPackedDepth);
        return;
    }

    vec2 offset1 = vec2(0.0, invTexSize.y);
    vec2 offset2 = vec2(invTexSize.x, 0.0);

    float selfDepth1 = getDepth(selfCoords + offset1);
    float selfDepth2 = getDepth(selfCoords + offset2);

    vec3 selfViewNormal = normalFromDepth(selfDepth, selfDepth1, selfDepth2, offset1, offset2);
    vec3 selfViewPos = screenSpaceToViewSpace(vec3(selfCoords, selfDepth), uInvProjection);

    vec3 randomVec = normalize(vec3(getNoiseVec2(selfCoords) * 2.0 - 1.0, 0.0));
    float pixelSize = getPixelSize(selfCoords, selfDepth);

    vec3 tangent = normalize(randomVec - selfViewNormal * dot(randomVec, selfViewNormal));
    vec3 bitangent = cross(selfViewNormal, tangent);
    mat3 TBN = mat3(tangent, bitangent, selfViewNormal);

    float occlusion = 0.0;
    for(int l = 0; l < dLevels; l++) {
        // TODO: smooth transition
        if (pixelSize * uDistanceFactor > uRadius[l]) continue;
        if (pixelSize * uMinDistanceFactor < uRadius[l]) continue;

        float levelOcclusion = 0.0;
        for(int i = 0; i < dNSamples; i++) {
            vec3 sampleViewPos = TBN * uSamples[i];
            sampleViewPos = selfViewPos + sampleViewPos * uRadius[l];

            vec4 offset = vec4(sampleViewPos, 1.0);
            offset = uProjection * offset;
            offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;

            float sampleDepth = getMappedDepth(offset.xy, selfCoords);
            if (uSolidBackground && sampleDepth == 1.0) {
                sampleViewPos = TBN * uSamples[i];
                sampleViewPos = selfViewPos + sampleViewPos * uRadius[l] * 0.5;

                offset = vec4(sampleViewPos, 1.0);
                offset = uProjection * offset;
                offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;
                sampleDepth = getMappedDepth(offset.xy, selfCoords);
            }
            float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, sampleDepth), uInvProjection).z;

            // float depth_difference = selfViewPos.z - sampleViewZ;
            // float rho = clamp((depth_difference - uRadius[l]) / depth_difference, 0.0, 1.0);
            // float rho = (depth_difference <= 0.0 || depth_difference > uRadius[l]) ? 1.0 : 0.0;
            // levelOcclusion += step(sampleViewPos.z + pixelSize, sampleViewZ) * rho * uBias[l];

            levelOcclusion += step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uRadius[l] / abs(selfViewPos.z - sampleViewZ)) * uBias[l];
        }
        occlusion = max(occlusion, levelOcclusion);
    }
    occlusion = 1.0 - (occlusion / float(dNSamples));

    vec2 packedOcclusion = packUnitIntervalToRG(clamp(occlusion, 0.01, 1.0));

    gl_FragColor = vec4(packedOcclusion, selfPackedDepth);
}
`;