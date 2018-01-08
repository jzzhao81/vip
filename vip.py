# -*- coding: utf-8 -*-
# @Author: jzzhao
# @Email:  jianzhou.zhao@gmail.com
# @Date:   2017-12-31 13:04:56
# @Last Modified by:   jzzhao
# @Last Modified time: 2018-01-07 23:20:45

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

orbit_dict = {
        's': ['s'],
        'p': ['py','pz','px'],
        'd': ['dxy', 'dyz', 'dz2', 'dxz', 'x2-y2'],
        'f': ['f-3','f-2','f-1','f0','f1','f2','f3'],
        't2g': ['dxy', 'dyz', 'dxz'],
        'eg': ['dz2', 'x2-y2']
        }

def rgbline(ax, kpt, eig, red, blue, green, alpha=1., lw=1.0):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb

    points = np.array([kpt, eig]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    nsegment = len(kpt) - 1
    r = np.abs( [0.5 * (red[i] + red[i + 1]) for i in range(nsegment)] )
    g = np.abs( [0.5 * (green[i] + green[i + 1]) for i in range(nsegment)] )
    b = np.abs( [0.5 * (blue[i] + blue[i + 1]) for i in range(nsegment)] )
    a = np.ones(nsegment, np.float) * alpha

    lc = LineCollection(segments, colors=list(zip(r, g, b, a)), linewidth=lw)
    ax.add_collection(lc)

class vasprun:

    def __init__(self, file=None):
        if file :
            self.xml_data = ET.parse(file)
        else :
            pass

    # get total dos
    # self.total_dos is a numpy.array, with [num_spin, num_data_grid, 3(energy, DOS, integrated DOS)]
    def get_dos(self, partial=False):
        self.dos = dos(xml_data=self.xml_data)
        self.dos.get_total_dos()
        if partial : self.dos.get_partial_dos()
        return

    def get_kpoints(self):
        self.kpoints = kpoints(xml_data=self.xml_data).get_klist()
        return

    def get_bandstructure(self,projected=False):
        self.bands = bandstructure(xml_data=self.xml_data)
        self.bands.get_eigenvalues()
        if projected : self.bands.get_eigenvalues_projected()
        return


class atom:

    def __init__(self, name='Atomic infomation'):
        self.name = name

    def from_xml(self, xml):

        # Load vasprun.xml file
        atominfo_tree = xml.find('.//atominfo')
        self.natom = int( atominfo_tree[0].text )
        self.ntype = int( atominfo_tree[1].text )
        self.elements = []
        self.elemass = []
        self.valence = []
        self.index = []
        self.orbit = []

        # get atom types infomation
        for itype in atominfo_tree.find('.//array[@name="atomtypes"]/set'):
            self.elements.append( itype[1].text.strip() )
            self.elemass.append( float(itype[2].text) )
            self.valence.append( float(itype[3].text) )
            self.index.append([])

        # get element index for position
        for ipos, ion in enumerate(atominfo_tree.find('.//array[@name="atoms"]/set')):
            itype = int( ion[1].text ) - 1
            self.index[itype].append(ipos)

        return self

class dos:

    def __init__(self, xml_data, name="Density of State"):
        self.xml_data = xml_data
        self.name = name
        self.nspin = int( self.xml_data.find('.//parameters/separator[@name="electronic"]/separator[@name="electronic spin"]/i[@name="ISPIN"]').text )

    def get_energy_grid(self):
        self.efermi = float( self.xml_data.find('.//dos/i[@name="efermi"]').text )
        self.energy_grid = np.array( [ data.text.split()[0] for data in  self.xml_data.find('.//dos/total/array/set')[0] ], dtype=np.float)
        return

    def get_total_dos(self):
        self.get_energy_grid()
        self.total_dos = np.array( [ [ data.text.split()[1] for data in ispin ] \
            for ispin in self.xml_data.find('.//dos/total/array/set') ], dtype=np.float)
        return

    def get_partial_dos(self):
        self.get_orbit_name()
        self.atoms = atom().from_xml(xml=self.xml_data)
        self.partial_dos = np.array( [ [ [ data.text.split()[1:] for data in spin] for spin in ion ] \
            for ion in self.xml_data.find('.//dos/partial/array/set') ], dtype=np.float )
        return

    def get_orbit_name(self):
        partial_tree = self.xml_data.findall('.//dos/partial/array/field')
        self.orbit_name = []
        for elem in partial_tree:
            self.orbit_name.append(elem.text.strip())
        if 'energy' in self.orbit_name : self.orbit_name.remove('energy')
        return

    def get_atom_index(self, key):
        return self.atoms.index[self.atoms.elements.index(key)]

    def get_orbit_index(self, key):
        return [self.orbit_name.index(item) for item in orbit_dict.get(key)]

    def get_plot(self, labl=None, xlim=None, ylim=None, lw=1.2, efermi=None, partial=None):

        if efermi : self.energy_grid -= efermi

        fig, ax = plt.subplots(figsize=(6.0, 5.0))

        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel("Density of States")

        # set x limit
        if not xlim :
            ax.set_xlim(self.energy_grid.min(), self.energy_grid.max())
        else :
            ax.set_xlim(xlim[0], xlim[1])

        # set y limit
        if not ylim :
            ax.set_ylim(self.total_dos.min(), self.total_dos.max())
        else :
            ax.set_ylim(ylim[0], ylim[1])

        ax.fill_between(self.energy_grid, 0.0, self.total_dos[0], color=(0.7, 0.7, 0.7),facecolor=(0.7, 0.7, 0.7))
        ax.plot( self.energy_grid, self.total_dos[0], lw=lw, label='Total DOS' )

        if partial :
            atoms_list = [ next(iter(atom.keys())) for atom in partial ]
            orbit_list = [ next(iter(atom.values())) for atom in partial ]
            for iatm, atom in enumerate(atoms_list):
                atom_index = self.get_atom_index(atom)
                for iorb, orbit in enumerate(orbit_list[iatm]):
                    orbit_index = self.get_orbit_index(orbit)
                    line = np.sum( np.sum( [self.partial_dos[iatm, 0, :, orbit_index] for iatm in atom_index], axis=0 ), axis=0 )
                    ax.plot( self.energy_grid, line ,lw=lw, \
                        label=next(iter(partial[iatm].keys()))+' '+ next(iter(partial[iatm].values()))[iorb])

        ax.legend()

        return plt


class kpoints:

    def __init__(self, xml_data, name='Vasp K-points'):
        self.xml_data = xml_data
        self.name = name

    def get_klist(self):

        # get reciprocal basis
        self.rec_basis = np.array( [ data.text.split() \
            for data in self.xml_data.find('.//structure[@name="finalpos"]/crystal/varray[@name="rec_basis"]')], dtype=np.float )
        # get K-Points mode
        self.kmode = self.xml_data.find('.//kpoints/generation').attrib["param"]
        if self.kmode == "listgenerated" :
            self.kdiv = int( self.xml_data.find('.//kpoints/generation/i[@name="divisions"]').text )
        else :
            self.kdiv = np.array( self.xml_data.find('.//kpoints/generation/v[@name="divisions"]').text.split(), dtype=np.int )
        self.klist = np.array( [ data.text.split() for data in self.xml_data.find('.//kpoints/varray[@name="kpointlist"]') ], dtype=np.float)
        self.ntotal = len(self.klist)
        self.nkseg = self.ntotal // self.kdiv
        return self

    def get_kpath(self):
        self.get_klist()
        self.kpath = np.array([0.0], dtype=np.float)
        for ikpt in range(1, len(self.klist)) :
            dist = np.dot((self.klist[ikpt] - self.klist[ikpt-1]),self.rec_basis)
            self.kpath = np.append(self.kpath, np.linalg.norm(dist)+self.kpath[-1])
        return self.kpath

class bandstructure:

    def __init__(self, xml_data, name='band structure object'):
        self.name = name
        self.xml_data = xml_data
        self.nbands = int( self.xml_data.find('.//parameters/separator[@name="electronic"]/i[@name="NBANDS"]').text )
        self.nspin = int( self.xml_data.find('.//parameters/separator[@name="electronic"]/separator[@name="electronic spin"]/i[@name="ISPIN"]').text )

    def get_eigenvalues(self):
        self.kpoints = kpoints(xml_data=self.xml_data)
        self.kpoints.get_kpath()
        self.efermi = float( self.xml_data.find('.//dos/i[@name="efermi"]').text )
        eigenvalues_tree = self.xml_data.find('.//calculation/eigenvalues/array/set')
        self.eigs = np.array( [ [ [ bnd.text.split()[0] for bnd in kpt ] for kpt in spin ] \
            for spin in eigenvalues_tree[:self.nspin] ], dtype=np.float ).reshape(-1,self.kpoints.ntotal,self.nbands)
        return

    def get_orbit_name(self):
        partial_tree = self.xml_data.findall('.//projected/array/field')
        self.orbit_name = []
        for elem in partial_tree:
            self.orbit_name.append(elem.text.strip())
        if 'energy' in self.orbit_name : self.orbit_name.remove('energy')
        return

    def get_eigenvalues_projected(self):
        self.get_orbit_name()
        self.atoms = atom().from_xml(xml=self.xml_data)
        proj_tree = self.xml_data.find('.//projected/array/set')
        self.proj = np.array( [ [[ [r.text.split() for r in bnd] for bnd in kpt] \
            for kpt in spin] for spin in proj_tree[:self.nspin] ], dtype=np.float )
        return

    def get_xtics(self):
        # get x ticks
        xtck = []
        xtck.extend( [self.kpoints.kpath[iksg*self.kpoints.kdiv] for iksg in range(self.kpoints.nkseg)]  )
        xtck.append( self.kpoints.kpath[-1] )
        return xtck

    def get_atom_index(self, key):
        return self.atoms.index[self.atoms.elements.index(key)]

    def get_orbit_index(self, key):
        return [self.orbit_name.index(item) for item in orbit_dict.get(key)]

    def get_normalize(self):
        norm = np.array( [ np.sum( self.proj[ispin,:,:,:,:], axis=(2,3) ) for ispin in range(self.nspin) ] )
        for ispin in range(self.nspin):
            for ikpt in range(self.kpoints.ntotal):
                for ibnd in range(self.nbands):
                    self.proj[ispin,ikpt,ibnd,:,:] /= norm[ispin,ikpt,ibnd]

    def get_plot(self, labl=None, xlim=None, ylim=None, lw=1.2, efermi=None, projected=None):

        if efermi : self.eigs -= efermi

        # get x ticks
        tics = self.get_xtics()

        fig, ax = plt.subplots(figsize=(6.0, 6.0))

        # set x limit
        if not xlim :
            ax.set_xlim(self.kpoints.kpath.min(), self.kpoints.kpath.max())
        else :
            ax.set_xlim(xlim[0], xlim[1])

        # set y limit
        if not ylim :
            ax.set_ylim(self.eigs.min(), self.eigs.max())
        else :
            ax.set_ylim(ylim[0], ylim[1])

        # plot band structure
        # First projection is red, Second projection is blue, others are green
        if projected :
            self.get_normalize()
            atoms_list = [ next(iter(atom.keys())) for atom in projected ]
            orbit_list = [ next(iter(atom.values())) for atom in projected ]
            num_proj = sum([len(item) for item in orbit_list])

            atom_index = [ self.get_atom_index(atom) for atom in atoms_list ]
            projection = []
            for ielem, elem in enumerate(atom_index):
                aaa = np.array( [ np.sum( self.proj[ispin,:,:,elem,:], axis=0) for ispin in range(self.nspin) ] )
                orbit_index = self.get_orbit_index( orbit_list[ielem][0] )
                bbb = np.array( [ np.sum( aaa[ispin,:,:,orbit_index], axis = 0 ) for ispin in range(self.nspin) ] )
                projection.append(bbb)
            rest = np.ones((self.nspin,self.kpoints.ntotal,self.nbands),dtype=np.float)  - projection[0] - projection[1]
            projection.append( rest )
            projection = np.array(projection)
            for ispin in range(self.nspin):
                for ibnd in range(self.nbands):
                    rgbline(ax, self.kpoints.kpath, self.eigs[ispin,:,ibnd], projection[0,ispin,:,ibnd], projection[1,ispin,:,ibnd], projection[2,ispin,:,ibnd])
        else :
            for ibnd in range(self.nbands):
                ax.plot(self.kpoints.kpath, self.eigs[0,:,ibnd], 'r-', linewidth=lw )
            if self.nspin == 2 :
                for ibnd in range(self.nbands):
                    ax.plot(self.kpoints.kpath, self.eigs[1,:,ibnd], 'b-', linewidth=lw )

        # Fermi level
        ax.text(self.kpoints.kpath.max(), 0.0, r'$E_F$')
        ax.plot([self.kpoints.kpath.min(), self.kpoints.kpath.max()], [0.0, 0.0],
             color='grey', linestyle='--', linewidth=0.8)

        # Set x ticks
        ax.set_xticks(tics)
        # Set x labels
        if labl:
            ax.set_xticklabels(labl)
        else:
            ax.set_xticklabels([])

        # plot high symmetry lines
        for itick in ax.get_xticks()[1:-1]:
            ax.plot([itick, itick], [ax.get_ylim()[0], ax.get_ylim()[1]],
                    linestyle='--', linewidth=0.8, color='grey')

        return ax


if __name__ == '__main__':

    # DOS
    ############################################################################
    # dosrun = vasprun(file='ggados.xml')
    # dosrun.get_dos(partial=True)
    # dosrun.dos.get_plot(xlim=[-6,4], ylim=[0, 20], efermi=dosrun.dos.efermi, \
    #     partial=[{'Zr':['d']},{'Te':['p']}] )
    # plt.show()
    # plt.savefig('ggados.png', dpi=300)

    # Bandstructure
    ############################################################################
    bsrun = vasprun(file='ggabnd.xml')
    bsrun.get_bandstructure(projected=True)
    bsrun.bands.get_plot(ylim=[-2, 2], efermi=bsrun.bands.efermi, \
        projected=[{'Te':'p'}, {'Zr':'d'}])
    plt.show()
