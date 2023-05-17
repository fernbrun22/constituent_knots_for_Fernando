from topoly import alexander
import numpy as np

class Protein():
    def __init__(self, filename, ndxes):
        self.name = filename[:4]
        self.coords = self.load_coords(filename)
        self.end_ndx = len(self.coords.keys())
        self.bridges = self.get_pairs(ndxes)
        self.constituent_knots = self.topol_loops()

    def load_coords(self, filename):
        coords = {}
        with open('input_data/{}'.format(filename), 'r') as f:
            for i, line in enumerate(f.readlines()):
                x,y,z = line.strip().split()
                coords[i+1] = (i+1, float(x), float(y), float(z))
        return coords

    def get_pairs(self, ndxes):
        pairs = []
        for i, ndx1 in enumerate(ndxes[:-1]):
            for j, ndx2 in enumerate(ndxes[i+1:]):
                pairs.append((ndx1,ndx2))
        return pairs

    def format_coords(self, coords):
        str_coords = []
        for i,x,y,z in coords:
            str_coords.append('{:d} {:.2f} {:.2f} {:.2f}'.format(i,x,y,z))
        return '\n'.join(str_coords)

    def topol_loops(self):
        constituent_knots = {}
        t0, t0_dir, t0_mass = self.topol_default_loop()
        with open('topology.txt', 'a+') as f:
            f.write('{}A\n'.format(self.name))
            f.write('Whole chain: probabilistic {}; direct {}; mass_center {}\n'.format(t0, t0_dir, t0_mass))
        constituent_knots['default'] = t0
        for bridge in self.bridges:
            print('\rPDB:{}, bridge:{:d}-{:d}'.format(self.name,bridge[0],bridge[1]),end='')
            t_in, dist = self.topol_internal_loop(bridge)
            with open('topology.txt', 'a+') as f:
                f.write('Bridge {:d}-{:d} dist: {:.2f}A\n'.format(bridge[0], bridge[1], dist))
                f.write('  Internal loop {:d}-{:d}: direct {}\n'.format(bridge[0], bridge[1], t_in))
            constituent_knots['in_{:d}-{:d}'.format(bridge[0],bridge[1])] = t_in
            t_sh, t_sh_dir, t_sh_mass = self.topol_shorter_loop(bridge)
            with open('topology.txt', 'a+') as f:
                f.write('  Shorter loop {:d}-{:d}: probabilistic {}; direct {}; mass_center {}\n'.format(bridge[0], bridge[1], t_sh, t_sh_dir, t_sh_mass))
            constituent_knots['sh_{:d}-{:d}'.format(bridge[0],bridge[1])] = t_sh
        with open('topology.txt', 'a+') as f:
            f.write('=====================\n')
        return constituent_knots
        
    def topol_default_loop(self):
        coords = [self.coords[i] for i in range(1,self.end_ndx+1)]
        with open('loop_xyzs/{}_default.xyz'.format(self.name), 'w') as f:
            f.write(self.format_coords(coords))
        return alexander(coords, max_cross=60, closure=2), alexander(coords, max_cross=60, closure=0), alexander(coords, max_cross=60, closure=1)

    def topol_internal_loop(self, bridge):
        coords = [self.coords[i] for i in range(bridge[0],bridge[1]+1)]
        with open('loop_xyzs/{}_in_{:d}-{:d}.xyz'.format(self.name,bridge[0],bridge[1]), 'w') as f:
            f.write(self.format_coords(coords))
        ndx1, x1, y1, z1 = self.coords[bridge[0]]
        ndx2, x2, y2, z2 = self.coords[bridge[1]+1]
        dist = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        return alexander(coords, max_cross=60, closure=0), dist

    def topol_shorter_loop(self, bridge):
        coords = [self.coords[i] for i in range(1,bridge[0]+1)]
        coords += [self.coords[i] for i in range(bridge[1],self.end_ndx+1)]
        with open('loop_xyzs/{}_sh_{:d}-{:d}.xyz'.format(self.name,bridge[0],bridge[1]), 'w') as f:
            f.write(self.format_coords(coords))
        return alexander(coords, max_cross=60, closure=2), alexander(coords, max_cross=60, closure=0), alexander(coords, max_cross=60, closure=1)


if __name__ == '__main__':
    auk = Protein('1aukA_CA.xyz', [12,263,264])
    ei6_1 = Protein('1ei6A_CA.xyz', [24,63,239,240,241])
    ei6_2 = Protein('1ei6A_CA.xyz', [201,205,367])
    ejj_1 = Protein('1ejjA_CA.xyz', [10,60,442,443])
    ejj_2 = Protein('1ejjA_CA.xyz', [401,405,460])

