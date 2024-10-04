/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { computeInteractions, Interactions, InteractionsParams as _InteractionsParams } from './interactions/interactions';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Features } from './interactions/features';
import { InteractionType, FeatureType } from './interactions/common';

export const InteractionsParams = {
    ..._InteractionsParams
};
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export type InteractionsValue = Interactions

export const InteractionsProvider: CustomStructureProperty.Provider<InteractionsParams, InteractionsValue> = CustomStructureProperty.createProvider({
    label: 'Interactions',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_interactions',
    }),
    type: 'local',
    defaultParams: InteractionsParams,
    getParams: (data: Structure) => InteractionsParams,
    isApplicable: (data: Structure) => !data.isCoarseGrained,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<InteractionsProps>) => {
        const p = { ...PD.getDefaultValues(InteractionsParams), ...props };
        const interactions = await computeInteractions(ctx, data, p);
        logRefinedInteractions(data, interactions);
        return { value: interactions };
    }
});

// Helper function to log data
function logData(filename: string, data: any) {
    const interactionData = JSON.stringify(data, getCircularReplacer(), 2);
    const blob = new Blob([interactionData], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

// Helper function to avoid circular references in JSON
function getCircularReplacer() {
    const seen = new WeakSet();
    return (key: any, value: any) => {
        if (typeof value === 'object' && value !== null) {
            if (seen.has(value)) {
                return;
            }
            seen.add(value);
        }
        return value;
    };
}

function logRefinedInteractions(structure: Structure, interactions: Interactions) {
    const dataToSave: any[] = [];
    const waterInteractionsMap = new Map<string, any[]>();

    const naturalAminoAcids = new Set([
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET',
        'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'
    ]);

    const getResidueName = (unit: Unit, index: number) => {
        return unit.model.atomicHierarchy.atoms.auth_comp_id.value(unit.elements[index]);
    };

    const getResidueID = (unit: Unit, index: number) => {
        const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[unit.elements[index]];
        return unit.model.atomicHierarchy.residues.auth_seq_id.value(residueIndex);
    };

    const getChainId = (unit: Unit, index: number) => {
        const chainIndex = unit.model.atomicHierarchy.chainAtomSegments.index[unit.elements[index]];
        return unit.model.atomicHierarchy.chains.label_asym_id.value(chainIndex);
    };

    const getAtomName = (unit: Unit, index: number) => {
        return unit.model.atomicHierarchy.atoms.label_atom_id.value(unit.elements[index]);
    };

    const getAromaticRingAtoms = (unit: Unit, featureIndex: number) => {
        const ringAtoms = [];
        const offsets = interactions.unitsFeatures.get(unit.id).offsets;
        const members = interactions.unitsFeatures.get(unit.id).members;

        for (let i = offsets[featureIndex]; i < offsets[featureIndex + 1]; i++) {
            const atomIndex = members[i];
            ringAtoms.push(getAtomName(unit, atomIndex));
        }
        return ringAtoms.join('-');
    };

    const isBackboneAtom = (atomName: string) => {
        return atomName === 'N' || atomName === 'C' || atomName === 'O' || atomName === 'CA';
    };

    const determineRole = (featureType: FeatureType) => {
        return featureType === FeatureType.HydrogenDonor || featureType === FeatureType.WeakHydrogenDonor ? 'Hydrogen Donor' : 'Hydrogen Acceptor';
    };

    const addInteraction = (key: string, interaction: any) => {
        if (!waterInteractionsMap.has(key)) {
            waterInteractionsMap.set(key, []);
        }

        // Check for duplicates before adding
        const existingInteractions = waterInteractionsMap.get(key)!;
        const isDuplicate = existingInteractions.some((existingInteraction) => {
            return (
                existingInteraction.unitA.residueName === interaction.unitA.residueName &&
                existingInteraction.unitA.residueID === interaction.unitA.residueID &&
                existingInteraction.unitB.residueName === interaction.unitB.residueName &&
                existingInteraction.unitB.residueID === interaction.unitB.residueID &&
                existingInteraction.unitA.role === interaction.unitA.role &&
                existingInteraction.unitB.role === interaction.unitB.role
            );
        });

        if (!isDuplicate) {
            waterInteractionsMap.get(key)!.push(interaction);
        }
    };

    for (let i = 0, il = interactions.contacts.edgeCount; i < il; ++i) {
        const e = interactions.contacts.edges[i];
        const uA = structure.unitMap.get(e.unitA) as Unit.Atomic;
        const uB = structure.unitMap.get(e.unitB) as Unit.Atomic;

        const infoA = Features.Info(structure, uA, interactions.unitsFeatures.get(e.unitA));
        infoA.feature = e.indexA;
        const infoB = Features.Info(structure, uB, interactions.unitsFeatures.get(e.unitB));
        infoB.feature = e.indexB;

        const elementA = infoA.members[infoA.offsets[infoA.feature]];
        const elementB = infoB.members[infoB.offsets[infoB.feature]];

        const residueNameA = getResidueName(uA, elementA);
        const residueNameB = getResidueName(uB, elementB);

        const interactionType = InteractionType[e.props.type];

        let atomInfoA = getAtomName(uA, elementA);
        let atomInfoB = getAtomName(uB, elementB);

        // Handle pi-stacking specifically
        if (interactionType === 'PiStacking') {
            atomInfoA = getAromaticRingAtoms(uA, e.indexA);
            atomInfoB = getAromaticRingAtoms(uB, e.indexB);
        }

        const residueA: any = {
            residueName: residueNameA,
            residueID: getResidueID(uA, elementA),
            chain: getChainId(uA, elementA),
            atom: atomInfoA,
        };

        const residueB: any = {
            residueName: residueNameB,
            residueID: getResidueID(uB, elementB),
            chain: getChainId(uB, elementB),
            atom: atomInfoB,
        };

        if (interactionType === 'HydrogenBond' || interactionType === 'WeakHydrogenBond') {
            // Only assign roles and locations if the interaction does not involve water
            if (residueNameA !== 'HOH' && residueNameB !== 'HOH') {
                residueA.role = determineRole(infoA.types[infoA.feature]);
                residueB.role = determineRole(infoB.types[infoB.feature]);

                if (naturalAminoAcids.has(residueNameA) && isBackboneAtom(residueA.atom)) {
                    residueA.location = 'Backbone';
                } else if (naturalAminoAcids.has(residueNameA)) {
                    residueA.location = 'Side chain';
                }

                if (naturalAminoAcids.has(residueNameB) && isBackboneAtom(residueB.atom)) {
                    residueB.location = 'Backbone';
                } else if (naturalAminoAcids.has(residueNameB)) {
                    residueB.location = 'Side chain';
                }
            }
        }

        const interaction = {
            interactionType,
            unitA: residueA,
            unitB: residueB
        };

        if (interactionType === 'HydrogenBond' || interactionType === 'WeakHydrogenBond') {
            if (residueNameA === 'HOH' || residueNameB === 'HOH') {
                const waterKey = residueNameA === 'HOH' ? `${residueNameA}-${residueA.residueID}-${residueA.chain}` : `${residueNameB}-${residueB.residueID}-${residueB.chain}`;
                addInteraction(waterKey, interaction);
            } else if (!naturalAminoAcids.has(residueNameA)) {
                dataToSave.push(interaction);
            }
        } else if (!naturalAminoAcids.has(residueNameA)) {
            dataToSave.push(interaction);
        }
    }

    // Grouping water-mediated hydrogen bonds
    waterInteractionsMap.forEach((interactions, key) => {
        if (interactions.length > 1) {
            const residues: any[] = [];
            let water: any = null;
            let ligand = null;

            interactions.forEach(interaction => {
                const otherUnit = interaction.unitA.residueName === 'HOH' ? interaction.unitB : interaction.unitA;
                if (interaction.unitA.residueName === 'HOH') {
                    water = interaction.unitA;
                } else if (interaction.unitB.residueName === 'HOH') {
                    water = interaction.unitB;
                }
                if (naturalAminoAcids.has(otherUnit.residueName)) {
                    residues.push(otherUnit);
                } else {
                    ligand = otherUnit;
                }
            });

            // Ensure no duplicate residues
            const uniqueResidues = Array.from(new Set(residues.map(res => JSON.stringify(res))))
                .map(res => JSON.parse(res));

            if (ligand && uniqueResidues.length > 0 && water) {
                const combinedInteraction = {
                    interactionType: 'Water-Bridged Hydrogen Bond',
                    residues: uniqueResidues,
                    water: {
                        residueName: 'HOH',
                        residueID: water.residueID,
                        chain: water.chain,
                        atom: water.atom
                    },
                    ligand: ligand
                };

                // Only add if it doesn't already exist
                if (!dataToSave.some(entry => JSON.stringify(entry) === JSON.stringify(combinedInteraction))) {
                    dataToSave.push(combinedInteraction);
                }
            }
        }
    });

    logData('refined_interactions.json', dataToSave);
}
