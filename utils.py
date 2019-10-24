import os
import re
import tempfile
from collections import OrderedDict, defaultdict
from copy import deepcopy

import freud
import gsd
import gsd.hoomd
import gsd.pygsd
import mbuild as mb
import numpy as np
from openbabel import pybel
from mbuild.exceptions import MBuildError
from mbuild.utils.io import import_
from oset import oset as OrderedSet
from parmed.periodic_table import Element


def bin_distribution(vals, nbins, start=None, stop=None):
    """
    Calculates a distribution given an array of data

    Parameters
    ----------
    vals : np.ndarry (N,), values over which to calculate the distribution
    start : float, value to start bins (default min(bonds_dict[bond]))
    stop : float, value to stop bins (default max(bonds_dict[bond]))
    step : float, step size between bins (default (stop-start)/30)

    Returns
    -------
    np.ndarray (nbins,2), where the first column is the mean value of the bin and
    the second column is number of values which fell into that bin
    """
    if start == None:
        start = min(vals)
    if stop == None:
        stop = max(vals)
    step = (stop-start)/nbins

    bins = [i for i in np.arange(start, stop, step)]
    dist = np.empty([len(bins)-1,2])
    for i, length in enumerate(bins[1:]):
        in_bin = [b for b in vals if b > bins[i] and b < bins[i+1]]
        dist[i,1] = len(in_bin)
        dist[i,0] = np.mean((bins[i],bins[i+1]))
    return dist


def autocorr1D(array):
    """
    Takes in a linear numpy array, performs autocorrelation
    function and returns normalized array with half the length
    of the input
    """
    ft = np.fft.rfft(array - np.average(array))
    acorr = np.fft.irfft(ft * np.conjugate(ft)) / (len(array) * np.var(array))
    return acorr[0 : len(acorr) // 2]


def get_decorr(acorr):
    """
    Returns the decorrelation time of the autocorrelation, a 1D numpy array
    """
    return np.argmin(acorr > 0)


def error_analysis(data):
    """
    Returns the standard and relative error given a dataset in a 1D numpy array
    """
    serr = np.std(data) / np.sqrt(len(data))
    rel_err = np.abs(100 * serr / np.average(data))
    return (serr, rel_err)


def get_angle(a, b, c):
    """
    Calculates the angle between three points a-b-c

    Parameters
    ----------
    a,b,c : np.ndarrays, positions of points a, b, and c

    Returns
    -------
    float, angle in radians
    """
    ba = a - b
    bc = c - b

    cos = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return np.arccos(cos)


def distance(pos1, pos2):
    """
    Calculates euclidean distance between two points.

    Parameters
    ----------
    pos1, pos2 : ((3,) numpy.ndarray), xyz coordinates
        (2D also works)

    Returns
    -------
    float distance
    """
    return np.linalg.norm(pos1 - pos2)


def v_distance(pos_array, pos2):
    """
    Calculates euclidean distances between all points in pos_array and pos2.

    Parameters
    ----------
    pos_array : ((N,3) numpy.ndarray), array of coordinates
    pos2 : ((3,) numpy.ndarray), xyz coordinate
        (2D also works)

    Returns
    -------
    (N,) numpy.ndarray of distances
    """
    return np.linalg.norm(pos_array - pos2, axis=1)


def get_molecules(snapshot):
    """
    Creates list of sets of connected atom indices

    This code adapted from Matias Thayer's:
    https://stackoverflow.com/questions/10301000/python-connected-components

    Parameters
    ----------
    snapshot : gsd.hoomd.Snapshot

    Returns
    -------
    list of sets of connected atom indices
    """

    def _snap_bond_graph(snapshot):
        """
        Given a snapshot from a trajectory create a graph of the bonds

        get_molecules(gsd.hoomd.Snapshot) --> dict of sets
        """
        bond_array = snapshot.bonds.group
        bond_graph = defaultdict(set)
        for row in bond_array:
            bond_graph[row[0]].add(row[1])
            bond_graph[row[1]].add(row[0])
        return bond_graph

    def _get_connected_group(node, already_seen):
        """
        This code adapted from Matias Thayer's:
        https://stackoverflow.com/questions/10301000/python-connected-components
        """

        result = set()
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_seen.add(node)
            nodes.update(graph[node] - already_seen)
            result.add(node)
        return result, already_seen

    graph = _snap_bond_graph(snapshot)

    already_seen = set()
    result = []
    for node in graph:
        if node not in already_seen:
            connected_group, already_seen = _get_connected_group(node, already_seen)
            result.append(connected_group)
    return result


def get_compound_rdf(compound, A_name, B_name, rmax=None, rdf=None):
    """
    This function calculates the radial distribution function given
    an mbuild compound, the names of the particles, and the dimensions of the box.

    Parameters
    ----------
    compound : CG_Compound
    A_name, B_name : str, name(s) of particle.name in compound
    rmax : float, maximum radius to consider. (default None)
        If none is given it'll be the minimum box length / 4
    rdf : freud.density.RDF, if provided, this function will accumulate an average rdf,
        otherwise it will provide the rdf only for the given compound. (default None)

    Returns
    -------
    freud.density.RDF
    """

    A_pos = compound.xyz[compound.get_name_inds(A_name), :]
    B_pos = compound.xyz[compound.get_name_inds(B_name), :]

    try:
        compound.box[0]
    except AttributeError(
        "No box found. Make sure you are using " + "CG_Compound and not mbuild.Compound"
    ):
        return
    except TypeError("Box has not been set"):
        return

    box = compound.box

    if rmax is None:
        rmax = min(box) / 4
    if rdf is None:
        rdf = freud.density.RDF(rmax=rmax, dr=0.01)

    box = freud.box.Box(*compound.box)

    rdf.accumulate(box, A_pos, B_pos)
    return rdf


def map_good_on_bad(good_mol, bad_mol):
    """
    This function takes a correctly-typed (good) and a poorly-typed (bad)
    pybel molecule and transfers the bond and atom typing from the good to
    the bad molecule but retains the atom positions.
    It assumes that both molecules have the same number of particles and
    they maintain their order.
    Changes:
    atom- type, isaromatic
    bond- order, isaromatic

    Parameters
    ----------
    good_mol, bad_mol : pybel.Molecule

    Returns
    -------
    pybel.Molecule
    """

    for i in range(1, good_mol.OBMol.NumAtoms()):
        good_atom = good_mol.OBMol.GetAtom(i)
        bad_atom = bad_mol.OBMol.GetAtom(i)
        bad_atom.SetType(good_atom.GetType())
        if good_atom.IsAromatic():
            bad_atom.SetAromatic()
        else:
            bad_atom.UnsetAromatic()

    for i in range(1, good_mol.OBMol.NumBonds()):
        good_bond = good_mol.OBMol.GetBond(i)
        bad_bond = bad_mol.OBMol.GetBond(i)
        bad_bond.SetBO(good_bond.GetBO())
        if good_bond.IsAromatic():
            bad_bond.SetAromatic()
        else:
            bad_bond.UnsetAromatic()

    return bad_mol


def save_mol_to_file(good_mol, filename):
    """
    This function takes a correctly-typed (good) pybel molecule and saves
    the bond and atom typing to a file for later use.

    Parameters
    ----------
    good_mol : pybel.Molecule
    filename : str, name of file

    use map_file_on_bad() to use this file
    """

    with open(filename, "w") as f:
        for i in range(1, good_mol.OBMol.NumAtoms()):
            good_atom = good_mol.OBMol.GetAtom(i)
            f.write(f"{good_atom.GetType()}   {good_atom.IsAromatic()}\n")

        for i in range(1, good_mol.OBMol.NumBonds()):
            good_bond = good_mol.OBMol.GetBond(i)
            f.write(f"{good_bond.GetBO()}   {good_bond.IsAromatic()}\n")


def map_file_on_bad(filename, bad_mol):
    """
    This function takes a filename containing correctly-typed and a poorly-typed (bad)
    pybel molecule and transfers the bond and atom typing from the good to
    the bad molecule but retains the atom positions.
    It assumes that both molecules have the same number of particles and
    they maintain their order.
    Changes:
    atom- type, isaromatic
    bond- order, isaromatic

    Parameters
    ----------
    filename : str, generated using save_mol_to_file()
    bad_mol : pybel.Molecule

    Returns
    -------
    pybel.Molecule
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    atoms = []
    bonds = []
    for line in lines:
        one, two = line.split()
        if not one.isdigit():
            atoms.append((one, two))
        else:
            bonds.append((one, two))

    for i in range(1, bad_mol.OBMol.NumAtoms()):
        bad_atom = bad_mol.OBMol.GetAtom(i)
        bad_atom.SetType(atoms[i - 1][0])
        if atoms[i - 1][1] == "True":
            bad_atom.SetAromatic()
        else:
            bad_atom.UnsetAromatic()

    for i in range(1, bad_mol.OBMol.NumBonds()):
        bad_bond = bad_mol.OBMol.GetBond(i)
        bad_bond.SetBO(int(bonds[i - 1][0]))
        if bonds[i - 1][1] == "True":
            bad_bond.SetAromatic()
        else:
            bad_bond.UnsetAromatic()

    return bad_mol


def has_number(string):
    """
    Returns True if string contains a number.
    Else returns False.
    """
    return bool(re.search("[0-9]", string))


def has_common_member(set_a, tup):
    """
    return True if set_a (set) and tup (tuple) share a common member
    else return False
    """
    set_b = set(tup)
    return set_a & set_b


def cg_comp(comp, bead_inds):
    """
    given an mbuild compound and bead_inds(list of tup)
    return coarse-grained mbuild compound
    """
    N = comp.n_particles

    for bead, smarts, bead_name in bead_inds:
        bead_xyz = comp.xyz[bead, :]
        avg_xyz = np.mean(bead_xyz, axis=0)
        bead = mb.Particle(name=bead_name, pos=avg_xyz)
        bead.smarts_string = smarts
        comp.add(bead)

    return comp, N


def cg_bonds(comp, beads, N):
    """
    add bonds based on bonding in aa compound
    return bonded mbuild compound
    """
    bonds = comp.get_bonds()
    bead_bonds = []
    for i, (bead_i, _, _) in enumerate(beads[:-1]):
        for j, (bead_j, _, _) in enumerate(beads[(i + 1) :]):
            for pair in bonds:
                if (pair[0] in bead_i) and (pair[1] in bead_j):
                    bead_bonds.append((i, j + i + 1))
                if (pair[1] in bead_i) and (pair[0] in bead_j):
                    bead_bonds.append((i, j + i + 1))
    for pair in bead_bonds:
        bond_pair = [
            particle
            for i, particle in enumerate(comp.particles())
            if i == pair[0] + N or i == pair[1] + N
        ]
        comp.add_bond(bond_pair)
    return comp


def num2str(num):
    """
    Returns a capital letter for positive integers up to 701
    e.g. num2str(0) = 'A'
    """
    if num < 26:
        return chr(num + 65)
    return "".join([chr(num // 26 + 64), chr(num % 26 + 65)])


def coarse(mol, bead_smarts, atomistic=False):
    """
    Creates a coarse-grained (CG) compound given a starting structure and
    smart strings for desired beads.

    Parameters
    ----------
    mol : pybel.Molecule
    bead_smarts : list of str, list of desired SMARTS strings of CG beads
    atomistic : bool, if True, the coarse-grain compound will be overlaid
    with the atomistic one (default=False)

    Returns
    -------
    CG_Compound
    """
    matches = []
    for i, smart_str in enumerate(bead_smarts):
        smarts = pybel.Smarts(smart_str)
        if not smarts.findall(mol):
            print(f"{smart_str} not found in compound!")
        for group in smarts.findall(mol):
            group = tuple(i - 1 for i in group)
            matches.append((group, smart_str, f"_{num2str(i)}"))

    seen = set()
    bead_inds = []
    for group, smarts, name in matches:
        # smart strings for rings can share atoms
        # add bead regardless of whether it was seen
        if has_number(smarts):
            for atom in group:
                seen.add(atom)
            bead_inds.append((group, smarts, name))
        # alkyl chains should be exclusive
        else:
            if has_common_member(seen, group):
                pass
            else:
                for atom in group:
                    seen.add(atom)
                bead_inds.append((group, smarts, name))

    n_atoms = mol.OBMol.NumHvyAtoms()
    if n_atoms != len(seen):
        print(
            "WARNING: Some atoms have been left out of coarse-graining!"
        )  # TODO make this more informative

    comp = CG_Compound.from_pybel(mol)
    comp, N = cg_comp(comp, bead_inds)
    comp = cg_bonds(comp, bead_inds, N)

    if not atomistic:
        comp.remove_atomistic()

    return comp


amber_dict = {
    "c": "C",
    "c1": "C",
    "c2": "C",
    "c3": "C",
    "ca": "C",
    "cp": "C",
    "cq": "C",
    "cc": "C",
    "cd": "C",
    "ce": "C",
    "cf": "C",
    "cg": "C",
    "ch": "C",
    "cx": "C",
    "cy": "C",
    "cu": "C",
    "cv": "C",
    "h1": "H",
    "h2": "H",
    "h3": "H",
    "h4": "H",
    "h5": "H",
    "ha": "H",
    "hc": "H",
    "hn": "H",
    "ho": "H",
    "hp": "H",
    "hs": "H",
    "hw": "H",
    "hx": "H",
    "f": "F",
    "cl": "Cl",
    "br": "Br",
    "i": "I",
    "n": "N",
    "n1": "N",
    "n2": "N",
    "n3": "N",
    "n4": "N",
    "na": "N",
    "nb": "N",
    "nc": "N",
    "nd": "N",
    "ne": "N",
    "nf": "N",
    "nh": "N",
    "no": "N",
    "o": "O",
    "oh": "O",
    "os": "O",
    "ow": "O",
    "p2": "P",
    "p3": "P",
    "p4": "P",
    "p5": "P",
    "pb": "P",
    "pc": "P",
    "pd": "P",
    "pe": "P",
    "pf": "P",
    "px": "P",
    "py": "P",
    "s": "S",
    "s2": "S",
    "s4": "S",
    "s6": "S",
    "sh": "S",
    "ss": "S",
    "sx": "S",
    "sy": "S",
}


# features SMARTS
features_dict = {
    "thiophene": "c1sccc1",
    "thiophene_F": "c1scc(F)c1",
    "alkyl_3": "CCC",
    "benzene": "c1ccccc1",
    "splitring1": "csc",
    "splitring2": "cc",
    "twobenzene": "c2ccc1ccccc1c2",
    "ring_F": "c1sc2c(scc2c1F)",
    "ring_3": "c3sc4cc5ccsc5cc4c3",
    "chain1": "OCC(CC)CCCC",
    "chain2": "CCCCC(CC)COC(=O)",
    "cyclopentadiene": "C1cccc1",
    "c4": "cC(c)(c)c",
    "cyclopentadienone": "C=C1C(=C)ccC1=O",
}


class CG_Compound(mb.Compound):
    def __init__(self):
        super().__init__()
        self.box = None

    @classmethod
    def from_gsd(cls, gsdfile, frame=-1, coords_only=False, scale=1.0):
        """
        Given a trajectory gsd file creates an mbuild.Compound.
        If there are multiple separate molecules, they are returned
        as one compound.

        Parameters
        ----------
        gsdfile : str, filename
        frame : int, frame number (default -1)
        coords_only : bool (default False)
            If True, return compound with no bonds
        scale : float, scaling factor multiplied to coordinates (default 1.0)

        Returns
        -------
        CG_Compound
        """
        f = gsd.pygsd.GSDFile(open(gsdfile, "rb"))
        t = gsd.hoomd.HOOMDTrajectory(f)

        snap = t[frame]
        bond_array = snap.bonds.group
        n_atoms = snap.particles.N

        # Add particles
        comp = cls()
        comp.box = snap.configuration.box[:3] * scale
        for i in range(n_atoms):
            name = snap.particles.types[snap.particles.typeid[i]]
            xyz = snap.particles.position[i] * scale
            charge = snap.particles.charge[i]

            atom = mb.Particle(name=name, pos=xyz, charge=charge)
            comp.add(atom, label=str(i))

        if not coords_only:
            # Add bonds
            for i in range(bond_array.shape[0]):
                atom1 = int(bond_array[i][0])
                atom2 = int(bond_array[i][1])
                comp.add_bond([comp[atom1], comp[atom2]])
        return comp

    def amber_to_element(self):
        """
        Pybel does not know how to parse atoms names in AMBER style
        so this functions renames them to their atomistic counterparts
        """
        for particle in self.particles():
            particle.name = amber_dict[particle.name]

    def remove_hydrogens(self):
        """
        Remove all particles with name = "H" in the compound
        """
        self.remove([i for i in self.particles() if i.name == "H"])

    def get_molecules(self):
        """
        Translates bond_graph.connected_components to particle indices in compound

        Returns
        -------
        list of sets of connected atom indices
        """
        particle_list = [part for part in self.particles()]
        molecules = []
        for group in self.bond_graph.connected_components():
            molecules.append(set(map(particle_list.index, group)))
        return molecules

    def get_bonds(self):
        """
        Translates bond_graph.bond_edges to particle indices in compound

        Returns
        -------
        list of tuples of bonded atom indices sorted
        """
        particle_list = [part for part in self.particles()]
        bonds = []
        for tup in self.bond_graph.edges():
            bonds.append(tuple(sorted(map(particle_list.index, tup))))
        # This sorting is required for coarse-graining
        bonds.sort(key=lambda tup: (tup[0], tup[1]))
        return bonds

    def from_pybel(pybel_mol, use_element=True):
        """
        Create a Compound from a Pybel.Molecule

        Parameters
        ---------
        pybel_mol: pybel.Molecule
        use_element : bool, default True
            If True, construct mb Particles based on the pybel Atom's element.
            If False, constructs mb Particles based on the pybel Atom's type

        Returns
        ------
        cmpd : CG_Compound
        """
        openbabel = mb.utils.io.import_("openbabel")
        cmpd = CG_Compound()
        resindex_to_cmpd = {}

        # Iterating through pybel_mol for atom/residue information
        # This could just as easily be implemented by
        # an OBMolAtomIter from the openbabel library,
        # but this seemed more convenient at time of writing
        # pybel atoms are 1-indexed, coordinates in Angstrom
        for atom in pybel_mol.atoms:
            xyz = np.array(atom.coords) / 10
            if use_element:
                try:
                    temp_name = Element[atom.atomicnum]
                except KeyError:
                    warn(
                        "No element detected for atom at index "
                        "{} with number {}, type {}".format(
                            atom.idx, atom.atomicnum, atom.type
                        )
                    )
                    temp_name = atom.type
            else:
                temp_name = atom.type
            temp = mb.compound.Particle(name=temp_name, pos=xyz)
            if hasattr(atom, "residue"):  # Is there a safer way to check for res?
                if atom.residue.idx not in resindex_to_cmpd:
                    res_cmpd = CG_Compound()
                    resindex_to_cmpd[atom.residue.idx] = res_cmpd
                    cmpd.add(res_cmpd)
                resindex_to_cmpd[atom.residue.idx].add(temp)
            else:
                cmpd.add(temp)

        # Iterating through pybel_mol.OBMol for bond information
        # Bonds are 0-indexed, but the atoms are 1-indexed
        # Bond information doesn't appear stored in pybel_mol,
        # so we need to look into the OBMol object,
        # using an iterator from the openbabel library
        for bond in openbabel.OBMolBondIter(pybel_mol.OBMol):
            cmpd.add_bond(
                [cmpd[bond.GetBeginAtomIdx() - 1], cmpd[bond.GetEndAtomIdx() - 1]]
            )

        if hasattr(pybel_mol, "unitcell"):
            box = mb.box.Box(
                lengths=[
                    pybel_mol.unitcell.GetA() / 10,
                    pybel_mol.unitcell.GetB() / 10,
                    pybel_mol.unitcell.GetC() / 10,
                ],
                angles=[
                    pybel_mol.unitcell.GetAlpha(),
                    pybel_mol.unitcell.GetBeta(),
                    pybel_mol.unitcell.GetGamma(),
                ],
            )
            cmpd.periodicity = box.lengths
        else:
            warn("No unitcell detected for pybel.Molecule {}".format(pybel_mol))
            box = None

        cmpd.box = box.maxs

        return cmpd

    def wrap(self):
        """
        Finds particles which are out of the box and translates
        them to within the box.
        """
        try:
            freud_box = freud.box.Box(*list(self.box))
        except TypeError:
            print("Can't wrap because CG_Compound.box values aren't assigned.")
            return
        particles = [part for part in self.particles()]
        # find rows where particles are out of the box
        for row in np.argwhere(abs(self.xyz) > self.box / 2)[:, 0]:
            new_xyz = freud_box.wrap(particles[row].pos)
            particles[row].translate_to(new_xyz)

    def unwrap(self, d_tolerance=0.22, _count=0):
        """
        Used to correct molecules which span the periodic boundary by translating particles
        to their real-space position. The function uses a distance tolerance to detect
        bonds which span the periodic boundary and from those determines which particle
        should be considered an outlier, finds other particles the outlier is bonded to,
        and shifts their position.

        Parameters
        ----------
        d_tolerance = float, distance beyond which a bond is considered "bad" (default=0.22)
        _count = int, used in recursive algorithm to prevent function from getting stuck
                 fixing bonds which span the pbc -- user should not set this value.

        if function is getting stuck in an endless loop, try adjusting d_tolerance
        """

        molecules = self.get_molecules()
        particles = [part for part in self.particles()]

        def check_bad_bonds(compound):
            """
            Used for identifying particles whose bonds span the periodic boundary.
            Finds particles indices in the compound with bonds longer than the
            distance tolerance.

            Parameters
            ----------
            compound : CG_Compound

            Returns
            -------
            list of tuples of particle indices
            """
            bad_bonds = [
                bond
                for bond in compound.bonds()
                if distance(bond[0].pos, bond[1].pos) > d_tolerance
            ]
            maybe_outliers = [
                (particles.index(bond[0]), particles.index(bond[1]))
                for bond in bad_bonds
            ]
            return maybe_outliers

        maybe_outliers = check_bad_bonds(self)
        if not maybe_outliers:
            print(
                f"No bonds found longer than {d_tolerance}. Either compound doesn't need"
                + " unwrapping or d_tolerance is too small. No changes made."
            )
            return

        def find_outliers(compound):
            """
            Finds "outliers" (bonded particles which span the periodic boundary).
            Starts by finding bonds that are too long, then determines which particle
            in the pair is an outlier based on whether removal of that particle reduces
            the average distance from the particles in the molecule to the geometric center.
            From these, the function follows the bond graph and adds all particles
            bonded to the outliers.

            Parameters
            ----------
            compound : CG_Compound

            Returns
            -------
            set of the particle indices of all outliers in the compound.
            """

            def _is_outlier(index):
                for molecule in molecules:
                    if index in molecule:
                        test_molecule = molecule.copy()
                        test_molecule.remove(index)
                        a = compound.xyz[list(molecule), :]
                        b = compound.xyz[list(test_molecule), :]
                        center_a = np.mean(a, axis=0)
                        center_b = np.mean(b, axis=0)
                        avg_dist_a = np.mean(v_distance(a, center_a))
                        avg_dist_b = np.mean(v_distance(b, center_b))
                return avg_dist_a > avg_dist_b

            def _d_to_center(index):
                for molecule in molecules:
                    if index in molecule:
                        mol_xyz = compound.xyz[list(molecule), :]
                        center = np.mean(mol_xyz, axis=0)
                        dist = distance(particles[index].pos, center)
                return dist

            outliers = set()
            checked = set()
            for tup in maybe_outliers:
                if _is_outlier(tup[0]) and _is_outlier(tup[1]):
                    # add whichever is further from the center
                    d_0 = _d_to_center(tup[0])
                    d_1 = _d_to_center(tup[1])
                    if d_0 > d_1:
                        outliers.add(tup[0])
                    elif d_1 > d_0:
                        outliers.add(tup[1])
                    else:
                        raise RuntimeError(
                            f"Can't determine which is outlier between indices {tup}"
                        )
                elif _is_outlier(tup[0]):
                    outliers.add(tup[0])
                elif _is_outlier(tup[1]):
                    outliers.add(tup[1])
                checked.add(tup[0])
                checked.add(tup[1])

            starts = outliers.copy()
            while starts:
                index = starts.pop()
                outliers.update(bond_dict[index] - checked)
                starts.update(bond_dict[index] - checked)
                checked.add(index)
            return outliers

        bond_dict = self.bond_dict()
        outliers = find_outliers(self)

        # organize outliers by molecule
        outlier_dict = defaultdict(set)
        for outlier in outliers:
            mol_ind = [i for i, mol in enumerate(molecules) if outlier in mol]
            if mol_ind:
                outlier_dict[mol_ind[0]].add(outlier)

        # find the center of the molecule without outliers
        for mol_ind in outlier_dict:
            molecule = molecules[mol_ind].copy()
            molecule -= outlier_dict[mol_ind]  # not outlier indices
            mol_xyz = self.xyz[list(molecule), :]
            mol_avg = np.mean(mol_xyz, axis=0)

            # translate the outlier to its real-space position found using
            # freud.box.unwrap. the direction is determined using the
            # difference between the particle position and the molecule center
            freud_box = freud.box.Box(*list(self.box))
            for outlier in outlier_dict[mol_ind]:
                image = mol_avg - particles[outlier].pos
                img = np.where(image > self.box / 2, 1, 0) + np.where(
                    image < -self.box / 2, -1, 0
                )
                new_xyz = freud_box.unwrap(particles[outlier].pos, img)
                particles[outlier].translate_to(new_xyz)

        # check if any bad bonds remain
        bad_bonds = check_bad_bonds(self)
        if bad_bonds:
            if _count < 5:
                _count += 1
                print(f"Bad bonds still present! Trying unwrap again. {_count}")
                self.unwrap(d_tolerance=d_tolerance, _count=_count)
            else:
                print("Bad bonds still present and try limit exceeded.")

    def bond_dict(self):
        """
        given an CG_Compound return an dict of the particle indices for each bond

        CG_Compound.bond_dict() --> dict of sets
        """
        parts = [part for part in self.particles()]
        bond_dict = defaultdict(set)
        for bond in self.bonds():
            bond_dict[int(parts.index(bond[0]))].add(int(parts.index(bond[1])))
            bond_dict[int(parts.index(bond[1]))].add(int(parts.index(bond[0])))
        return bond_dict

    def get_name_inds(self, name):
        """
        Find indices of particles in compound where particle.name matches given name

        Parameters
        ----------
        name : str, particle.name in mb.Compound

        Returns
        -------
        list of particles indices which match name
        """
        return [i for i, part in enumerate(self.particles()) if part.name == name]

    def tuple_to_names(self, tup):
        """
        Get the names of particle indices passed in as a tuple.

        Parameters
        ----------
        tup : tuple of ints, particle indices

        Returns
        -------
        tuple of strings, particle.name of given indices
        """
        particles = [part for part in self.particles()]

        types = []
        for index in tup:
            types.append(particles[index].name)
        return tuple(sorted(types))

    def find_angles(self):
        """
        Adapted from cme_utils.manip.builder.building_block.find_angles()
        Finds unique angle constraints and their types.

        Returns
        -------
        Dictionary with keys correponding to the angle types and
        values which list the particle indices which have this angle type
        """
        angles = []
        bond_dict = self.bond_dict()
        for i in range(self.n_particles):
            for n1 in bond_dict[i]:
                for n2 in bond_dict[n1]:
                    if n2 != i:
                        if n2 > i:
                            angles.append((i, n1, n2))
                        else:
                            angles.append((n2, n1, i))
        angles = sorted(set(angles))
        angle_types = []
        for t in angles:
            angle_types.append(self.tuple_to_names(t))

        angle_dict = defaultdict(list)
        for a, b in zip(angle_types, angles):
            angle_dict[a].append(b)
        return angle_dict

    def find_bonds(self):
        """
        Finds unique bond constraints and their types.

        Returns
        -------
        Dictionary with keys correponding to the bond types and
        values which list the particle indices which have this type
        """
        bonds = []
        bond_dict = self.bond_dict()
        for i in range(self.n_particles):
            for n1 in bond_dict[i]:
                if n1 > i:
                    bonds.append((i, n1))
                else:
                    bonds.append((n1, i))
        bonds = sorted(set(bonds))
        bond_types = []
        for t in bonds:
            bond_types.append(self.tuple_to_names(t))

        bond_dict = defaultdict(list)
        for a, b in zip(bond_types, bonds):
            bond_dict[a].append(b)
        return bond_dict

    def find_pairs(self):
        """
        Finds unique (coarse-grained) pair types
        (coarse particle names start with "_")

        Returns
        -------
        list of tuples of pair names
        """
        particles = {p.name for p in self.particles() if p.name[0] == "_"}
        pairs = set()
        for i in particles:
            for j in particles:
                pair = tuple(sorted([i, j]))
                pairs.add(pair)
        return sorted(pairs)

    def remove_atomistic(self):
        """
        all coarse-grained particles are named starting with '_'
        remove any particles whose names do not start with '_'
        """
        for i in self.particles():
            if i.name[0] != "_":
                self.remove(i)
        # Remove residual ports
        for child in self.children:
            if type(child) == mb.port.Port:
                self.remove(child)

    def is_bad_bond(self, tup):
        """
        Determines whether a bond spans the periodic boundary based on a distance
        cutoff of the self.box/2

        Parameters
        ----------
        tup : tuple, indices of the bonded particles

        Returns
        -------
        bool
        """
        if tup not in self.get_bonds() and tup[::-1] not in self.get_bonds():
            print(f"Bond {tup} not found in compound! Aborting...")
            return
        pair = [p for i, p in enumerate(self.particles()) if i == tup[0] or i == tup[1]]
        test = np.where(abs(pair[0].xyz - pair[1].xyz) > self.box / 2)[1]
        if test.size > 0:
            return True
        else:
            return False

    def unwrap_position(self, tup):
        """
        Given the indices of a bonded pair which spans the periodic boundary,
        moves the second index to it's real-space position.

        Parameters
        ----------
        tup : tuple, indices (2) of bonded particles

        Returns
        -------
        np.ndarray(3,), unwrapped coordinates for index in tup[1]
        (if you want to move the first index, enter it as tup[::-1])
        """
        freud_box = freud.box.Box(*list(self.box))
        pair = [p for i, p in enumerate(self.particles()) if i == tup[0] or i == tup[1]]
        diff = pair[0].pos - pair[1].pos
        img = np.where(diff > self.box[:3] / 2, 1, 0) + np.where(
            diff < -self.box[:3] / 2, -1, 0
        )
        return freud_box.unwrap(pair[1].pos, img)


    @classmethod
    def from_mbuild(cls, compound):
        """
        Instantiates a CG_Compound and follows mb.Compound.deep_copy
        to copy particles and bonds to CG_Compound

        Parameters
        ----------
        compound : mb.Compound to be compied

        Returns
        -------
        CG_Compound
        """

        comp = cls()

        clone_dict = {}
        comp.name = deepcopy(compound.name)
        comp.periodicity = deepcopy(compound.periodicity)
        comp._pos = deepcopy(compound._pos)
        comp.port_particle = deepcopy(compound.port_particle)
        comp._check_if_contains_rigid_bodies = deepcopy(
            compound._check_if_contains_rigid_bodies
        )
        comp._contains_rigid = deepcopy(compound._contains_rigid)
        comp._rigid_id = deepcopy(compound._rigid_id)
        comp._charge = deepcopy(compound._charge)

        if compound.children is None:
            comp.children = None
        else:
            comp.children = OrderedSet()
        # Parent should be None initially.
        comp.parent = None
        comp.labels = OrderedDict()
        comp.referrers = set()
        comp.bond_graph = None
        for p in compound.particles():
            new_particle = mb.Particle(name=p.name, pos=p.xyz.flatten())
            comp.add(new_particle)
            clone_dict[p] = new_particle

        for c1, c2 in compound.bonds():
            try:
                comp.add_bond((clone_dict[c1], clone_dict[c2]))
            except KeyError:
                raise MBuildError(
                    "Cloning failed. Compound contains bonds to "
                    "Particles outside of its containment hierarchy."
                )
        return comp


    def _visualize_py3dmol(self, show_ports=False, color_scheme={}):
        """
        Visualize the Compound using py3Dmol.
        Allows for visualization of a Compound within a Jupyter Notebook.
        Modified to show atomistic elements (translucent) with larger CG beads.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}
        Returns
        ------
        view : py3Dmol.view
        """
        py3Dmol = import_("py3Dmol")

        coarse = mb.clone(self)
        atomistic = mb.clone(self)
        atomistic.remove_coarse()
        coarse.remove_atomistic()

        modified_color_scheme = {}
        for name, color in color_scheme.items():
            # Py3dmol does some element string conversions,
            # first character is as-is, rest of the characters are lowercase
            new_name = name[0] + name[1:].lower()
            modified_color_scheme[new_name] = color
            modified_color_scheme[name] = color

        cg_names = []
        for particle in coarse.particles():
            if not particle.name:
                particle.name = "UNK"
            else:
                if (particle.name != 'Compound') and (particle.name != 'CG_Compound'):
                    cg_names.append(particle.name)

        atom_names = []
        for particle in atomistic.particles():
            if not particle.name:
                particle.name = "UNK"
            else:
                if (particle.name != 'Compound') and (particle.name != 'CG_Compound'):
                    atom_names.append(particle.name)

        tmp_dir = tempfile.mkdtemp()

        view = py3Dmol.view()

        if atom_names:
            atomistic.save(
                os.path.join(tmp_dir, "atomistic_tmp.mol2"),
                show_ports=show_ports,
                overwrite=True,
            )

            # atomistic
            with open(os.path.join(tmp_dir, "atomistic_tmp.mol2"), "r") as f:
                view.addModel(f.read(), "mol2", keepH=True)

            if cg_names:
                opacity = 0.6
            else:
                opacity = 1.0

            view.setStyle(
                {
                    "stick": {"radius": 0.2, "opacity": opacity, "color": "grey"},
                    "sphere": {
                        "scale": 0.3,
                        "opacity": opacity,
                        "colorscheme": modified_color_scheme,
                    },
                }
            )

        # coarse
        if cg_names:
            coarse.save(
                os.path.join(tmp_dir, "coarse_tmp.mol2"),
                show_ports=show_ports,
                overwrite=True,
            )
            with open(os.path.join(tmp_dir, "coarse_tmp.mol2"), "r") as f:
                view.addModel(f.read(), "mol2", keepH=True)

            view.setStyle(
                {"atom": cg_names},
                {
                    "stick": {"radius": 0.2, "opacity": 1, "color": "grey"},
                    "sphere": {
                        "scale": 0.7,
                        "opacity": 1,
                        "colorscheme": modified_color_scheme,
                    },
                },
            )

        view.zoomTo()

        return view

    def remove_coarse(self):
        """
        all coarse-grained particles are named starting with '_'
        remove any particles whose names start with '_'
        """
        for i in self.particles():
            if i.name[0] == "_":
                self.remove(i)
        # Remove residual ports
        for child in self.children:
            if type(child) == mb.port.Port:
                self.remove(child)
