/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Quat, Vec3 } from '../mol-math/linear-algebra';
import { degToRad } from '../mol-math/misc';
import { Loci } from '../mol-model/loci';
import { StructureElement, Structure } from '../mol-model/structure';
import { Script } from '../mol-script/script';
import { BitFlags } from '../mol-util/bit-flags';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { stringToWords } from '../mol-util/string';

export { Clipping };

type Clipping = {
    readonly layers: ReadonlyArray<Clipping.Layer>
    readonly variant: Clipping.Variant
    readonly objects: Clipping.Objects
}

function Clipping(layers: Clipping['layers'], variant: Clipping['variant'], objects: Clipping['objects']): Clipping {
    return { layers, variant, objects };
}

function createClipObjects() {
    return {
        count: 0,
        type: (new Array(5)).fill(1),
        invert: (new Array(5)).fill(false),
        position: (new Array(5 * 3)).fill(0),
        rotation: (new Array(5 * 4)).fill(0),
        scale: (new Array(5 * 3)).fill(1),
    };
}

namespace Clipping {
    export type Layer = { readonly loci: StructureElement.Loci, readonly groups: Groups }
    export const Empty: Clipping = {
        layers: [],
        variant: 'pixel',
        objects: createClipObjects()
    };

    export type Groups = BitFlags<Groups.Flag>
    export namespace Groups {
        export const is: (g: Groups, f: Flag) => boolean = BitFlags.has;
        export const enum Flag {
            None = 0x0,
            One = 0x1,
            Two = 0x2,
            Three = 0x4,
            Four = 0x8,
            Five = 0x10,
            Six = 0x20,
        }

        export function create(flags: Flag): Groups {
            return BitFlags.create(flags);
        }

        export const Names = {
            'one': Flag.One,
            'two': Flag.Two,
            'three': Flag.Three,
            'four': Flag.Four,
            'five': Flag.Five,
            'six': Flag.Six,
        };
        export type Names = keyof typeof Names

        export function isName(name: string): name is Names {
            return name in Names;
        }

        export function fromName(name: Names): Flag {
            switch (name) {
                case 'one': return Flag.One;
                case 'two': return Flag.Two;
                case 'three': return Flag.Three;
                case 'four': return Flag.Four;
                case 'five': return Flag.Five;
                case 'six': return Flag.Six;
            }
        }

        export function fromNames(names: Names[]): Flag {
            let f = Flag.None;
            for (let i = 0, il = names.length; i < il; ++i) {
                f |= fromName(names[i]);
            }
            return f;
        }

        export function toNames(groups: Groups): Names[] {
            const names: Names[] = [];
            if (is(groups, Flag.One)) names.push('one');
            if (is(groups, Flag.Two)) names.push('two');
            if (is(groups, Flag.Three)) names.push('three');
            if (is(groups, Flag.Four)) names.push('four');
            if (is(groups, Flag.Five)) names.push('five');
            if (is(groups, Flag.Six)) names.push('six');
            return names;
        }
    }

    /** Clip object types */
    export const Type = {
        none: 0, // to switch clipping off
        plane: 1,
        sphere: 2,
        cube: 3,
        cylinder: 4,
        infiniteCone: 5,
    };

    export type Variant = 'instance' | 'pixel'

    export type Objects = {
        count: number
        type: number[]
        invert: boolean[]
        position: number[]
        rotation: number[]
        scale: number[]
    }

    export const Params = {
        variant: PD.Select('instance', PD.arrayToOptions<Variant>(['instance', 'pixel'])),
        objects: PD.ObjectList({
            type: PD.Select('plane', PD.objectToOptions(Type, t => stringToWords(t))),
            invert: PD.Boolean(false),
            position: PD.Vec3(Vec3()),
            rotation: PD.Group({
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { description: 'Angle in Degrees' }),
            }, { isExpanded: true }),
            scale: PD.Vec3(Vec3.create(1, 1, 1)),
        }, o => stringToWords(o.type))
    };
    export type Params = typeof Params
    export type Props = PD.Values<Params>

    export type Clip = {
        variant: Clipping['variant'],
        objects: Clipping['objects']
    }

    const qA = Quat();
    const qB = Quat();
    const vA = Vec3();
    const vB = Vec3();

    export function getClip(props: Props, clip?: Clip): Clip {
        const { type, invert, position, rotation, scale } = clip?.objects || createClipObjects();
        for (let i = 0, il = props.objects.length; i < il; ++i) {
            const p = props.objects[i];
            type[i] = Type[p.type];
            invert[i] = p.invert;
            Vec3.toArray(p.position, position, i * 3);
            Quat.toArray(Quat.setAxisAngle(qA, p.rotation.axis, degToRad(p.rotation.angle)), rotation, i * 4);
            Vec3.toArray(p.scale, scale, i * 3);
        }
        return {
            variant: props.variant,
            objects: { count: props.objects.length, type, invert, position, rotation, scale }
        };
    }

    export function areEqual(cA: Clipping, cB: Clipping) {
        if (cA.layers.length !== cB.layers.length) return false;
        for (let i = 0, il = cA.layers.length; i < il; ++i) {
            if (cA.layers[i].groups !== cB.layers[i].groups) return false;
            if (!Loci.areEqual(cA.layers[i].loci, cB.layers[i].loci)) return false;
        }
        if (cA.variant !== cB.variant) return false;
        if (cA.objects.count !== cB.objects.count) return false;

        const oA = cA.objects, oB = cB.objects;
        for (let i = 0, il = oA.count; i < il; ++i) {
            if (oA.invert[i] !== oB.invert[i]) return false;
            if (oA.type[i] !== oB.type[i]) return false;

            Vec3.fromArray(vA, oA.position, i * 3);
            Vec3.fromArray(vB, oB.position, i * 3);
            if (!Vec3.equals(vA, vB)) return false;

            Vec3.fromArray(vA, oA.scale, i * 3);
            Vec3.fromArray(vB, oB.scale, i * 3);
            if (!Vec3.equals(vA, vB)) return false;

            Quat.fromArray(qA, oA.rotation, i * 4);
            Quat.fromArray(qB, oB.rotation, i * 4);
            if (!Quat.equals(qA, qB)) return false;
        }
        return true;
    }

    /** Check if layers empty */
    export function isEmpty(clipping: Clipping) {
        return clipping.layers.length === 0;
    }

    /** Remap layers */
    export function remap(clipping: Clipping, structure: Structure): Clipping {
        const layers: Clipping.Layer[] = [];
        for (const layer of clipping.layers) {
            let { loci, groups } = layer;
            loci = StructureElement.Loci.remap(loci, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers, variant: clipping.variant, objects: clipping.objects };
    }

    /** Merge layers */
    export function merge(clipping: Clipping): Clipping {
        if (isEmpty(clipping)) return clipping;
        const { structure } = clipping.layers[0].loci;
        const map = new Map<Groups, StructureElement.Loci>();
        let shadowed = StructureElement.Loci.none(structure);
        for (let i = 0, il = clipping.layers.length; i < il; ++i) {
            let { loci, groups } = clipping.layers[il - i - 1]; // process from end
            loci = StructureElement.Loci.subtract(loci, shadowed);
            shadowed = StructureElement.Loci.union(loci, shadowed);
            if (!StructureElement.Loci.isEmpty(loci)) {
                if (map.has(groups)) {
                    loci = StructureElement.Loci.union(loci, map.get(groups)!);
                }
                map.set(groups, loci);
            }
        }
        const layers: Clipping.Layer[] = [];
        map.forEach((loci, groups) => {
            layers.push({ loci, groups });
        });
        return { layers, variant: clipping.variant, objects: clipping.objects };
    }

    /** Filter layers */
    export function filter(clipping: Clipping, filter: Structure): Clipping {
        if (isEmpty(clipping)) return clipping;
        const { structure } = clipping.layers[0].loci;
        const layers: Clipping.Layer[] = [];
        for (const layer of clipping.layers) {
            let { loci, groups } = layer;
            // filter by first map to the `filter` structure and
            // then map back to the original structure of the clipping loci
            const filtered = StructureElement.Loci.remap(loci, filter);
            loci = StructureElement.Loci.remap(filtered, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers, variant: clipping.variant, objects: clipping.objects };
    }

    export type ScriptLayer = { script: Script, groups: Groups }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure, clip: Clip): Clipping {
        const layers: Clipping.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, groups } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers, variant: clip.variant, objects: clip.objects };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, groups: Groups }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure, clip: Clip): Clipping {
        const layers: Clipping.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, groups } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, groups });
        }
        return { layers, variant: clip.variant, objects: clip.objects };
    }

    export function toBundle(clipping: Clipping) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = clipping.layers.length; i < il; ++i) {
            const { loci, groups } = clipping.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, groups });
        }
        return { layers };
    }
}