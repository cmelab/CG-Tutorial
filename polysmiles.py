import deepsmiles
import mbuild as mb


def convert_smiles(smiles=False, deep=False):
    """
    smiles and deep must be str format
    Converts from SMILES to DeepSMILES and vice versa.
    Whichever has a string provided, will convert to the other.
    If strings are proivded for both, then nothing happens
    """

    converter = deepsmiles.Converter(rings=True, branches=True)
    if smiles and deep:
        print("Only provide a string for one of smiles or deep")
        return ()
    if smiles:  # Convert from SMILES to DeepSMILES
        deep_string = converter.encode(smiles)
        return deep_string
    if deep:  # Convert from DeepSMILES to SMILES
        smiles_string = converter.decode(deep)
        return smiles_string


def poly_smiles(monomer_string, length=2):
    """
    Builds a polymer of desired length when given the deepsmiles string of the monomer.
    The smiles string must include asterisks, **, around the desired polymer site
      e.g., ccc*c*cc6

    Parameters
    ----------
    monomer_string : str, deepsmiles format of repeated structure
    length : int, number of times to repeat the monomer
    """

    # Find how many brackets are required at polymerization site
    atom_count = 0
    bracket_count = 0
    for s in monomer_string:
        if s.isalpha():
            atom_count += 1
        if s == ")":
            bracket_count += 1
    if bracket_count == 0:
        brackets = ")" * atom_count
    elif bracket_count != 0:
        brackets = ")" * (atom_count - bracket_count)

    # Find index num of poly site on modified DEEP SMILES string
    monomer_list = list(monomer_string)
    if "*" not in monomer_list:
        return print("ERROR: Identify the wanted polymerization site using *x*")
    key_indices = [i for i, value in enumerate(monomer_list) if value == "*"]

    # Checks for only a single given poly site
    if len(key_indices) != 2:
        return print("ERROR: Select only one polymerization site using *x*")

    # Check that the * are surrounding only a single atom
    if (key_indices[1] - key_indices[0] != 2):
        return print("ERROR: Select only one polymerization site using *x*")

    # Create poly site+brackets to the right of the atom
    monomer_list[key_indices[1]] = "{}" + "{}".format(brackets)
    monomer_list.remove("*")

    # Monomer string with the needed {} and without second * in the string
    template = "".join(monomer_list)
    monomer_list.remove("{}" + "{}".format(brackets))

    # Pure deepsmiles monomer string without {} or *
    monomer = "".join(monomer_list)
    # Loop & format polymer
    polymer = "{}"
    for i in range(length):
        if i == length - 1:
            polymer = polymer.format(monomer)
            break
        polymer = polymer.format(template)
    polymer_smiles = convert_smiles(deep=polymer)

    return polymer_smiles
