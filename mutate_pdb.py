import pdbfixer
import openmm.app as app
from openmm.app import PDBFile
import tempfile
def pdb_mutate(path_data, pdb_id, chain, mut_info, mutseq_pos, mut_resid):
    """This function mutates the wild residue type to mutant residue and fix the pdb block accordingly. This function
    uses PDBFixer """
    temp_file = path_data + pdb_id + '.af2.pdb'
    fixer0 = pdbfixer.PDBFixer(filename=temp_file)

    fixer0.findMissingResidues()
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer0.topology.chains())
    keys = fixer0.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chainA = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chainA.residues()))):
            missingResidues[key] = fixer0.missingResidues[key]
    fixer0.missingResidues = missingResidues

    fixer0.findMissingAtoms()
    fixer0.addMissingAtoms()
    PDBFile.writeFile(fixer0.topology, fixer0.positions, open(temp_file, 'w'), keepIds=True)
    #print(pdb_id.lower(), chain, mut_info)
    temp_mutfile = path_data + pdb_id + '_' + str(mutseq_pos) + mut_resid+'.af2.pdb'
    #print(temp_mutfile)
    '''if os.path.exists(temp_mutfile):
        return'''
    fixer = pdbfixer.PDBFixer(filename=path_data + pdb_id.lower() + '.af2.pdb')
    fixer.applyMutations([mut_info], chain)

    fixer.findMissingResidues()
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missingResidues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missingResidues

    #fixer.findMissingAtoms()
    #fixer.addMissingAtoms()
    with tempfile.NamedTemporaryFile(mode='w+') as temp_pdb:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, temp_pdb)
        temp_pdb.flush()
    PDBFile.writeFile(fixer.topology, fixer.positions, open(temp_mutfile, 'w'), keepIds=True)
    '''mol = Chem.MolFromPDBFile(temp_mutfile)
    mol = pdbmol_prepare.PreparePDBMol(mol, remove_incomplete=True, add_missing_atoms=True, removeHs=False,
                                       removeHOHs=False, residue_blacklist=[translate_1aa3(mut_resid)])
    #Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    #Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    mol = mol_fix.fix_mol(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    Chem.MolToPDBFile(mol, temp_mutfile)'''