from topoly import alexander
import numpy as np

class Protein():
    def __init__(self, filename, ion_neighbours, ss_bonds=[]):
        self.name = filename[:4]
        self.coords = self.load_coords(filename)
        self.end_ndx = len(self.coords.keys())
        self.ss_bridges = ss_bonds
        self.ion_bridges = self.get_pairs(ion_neighbours)
        self.bridges = sorted(self.ss_bridges + self.ion_bridges)
        with open('topology.txt', 'a+') as f:
            f.write('======================================\n')
            f.write('{}A\nSS bridges: {}\nION bridges: {}\n'.format(self.name,str(self.ss_bridges),str(self.ion_bridges)))
        self.topol_loops_single()
        self.topol_loops_double()

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

        for x,y,z in coords:
            str_coords.append('{:.2f} {:.2f} {:.2f}'.format(x,y,z))
        return '\n'.join(str_coords)

    def topol_loops_single(self):
        constituent_knots = {}
        aaaa = self.topol_default_loop()
        print(aaaa)
        t0, t0_dir, t0_mass = aaaa
        with open('topology.txt', 'a+') as f:
            f.write('Whole chain: probabilistic {}; direct {}; mass_center {}\n'.format(t0, t0_dir, t0_mass))
            f.write('=====================\n')
        constituent_knots['default'] = t0
        for bridge in self.bridges:
            print('\rPDB:{}, bridge:{:d}-{:d}'.format(self.name,bridge[0],bridge[1]),end='')
            t_in = self.topol_internal_loop(bridge)
            with open('topology.txt', 'a+') as f:
                f.write('Bridge {:d}-{:d}\n'.format(bridge[0], bridge[1]))
                f.write('  Internal loop {:d}-{:d}: direct {}\n'.format(bridge[0], bridge[1], t_in))
            constituent_knots['in_{:d}-{:d}'.format(bridge[0],bridge[1])] = t_in
            t_sh, t_sh_dir, t_sh_mass = self.topol_shorter_loop(bridge)
            with open('topology.txt', 'a+') as f:
                f.write('  Shorter loop {:d}-{:d}: probabilistic {}; direct {}; mass_center {}\n'.format(bridge[0], bridge[1], t_sh, t_sh_dir, t_sh_mass))
            constituent_knots['sh_{:d}-{:d}'.format(bridge[0],bridge[1])] = t_sh
            with open('topology.txt', 'a+') as f:
                f.write('=====================\n')
        print(constituent_knots)
        return constituent_knots

    def topol_loops_double(self):
        print(self.bridges)
        constituent_knots = {}
        for i,bridge1 in enumerate(self.bridges[:-1]):
            for j,bridge2 in enumerate(self.bridges[i+1:]):
                if len(set([*bridge1,*bridge2])) < 4:
                    # degeneration
                    topol_type = 'DEGENERATION'
                elif bridge1[1] < bridge2[0]:
                    # handcuff
                    topol_type = 'HANDCUFF'
                elif bridge1[1] > bridge2[1]:
                    # snowman
                    topol_type = 'SNOWMAN'
                elif bridge1[1] > bridge2[0]:
                    # theta
                    topol_type = 'THETA'
                else:
                    raise
                print('\rPDB:{}, {}, bridges:{:d}-{:d},{:d}-{:d}'.format(self.name,topol_type,bridge1[0],bridge1[1],bridge2[0],bridge2[1]),end='')
                with open('topology.txt', 'a+') as f:
                    f.write('{}, {}, bridges:{:d}-{:d},{:d}-{:d}\n'.format(self.name,topol_type,bridge1[0],bridge1[1],bridge2[0],bridge2[1]))
                if topol_type == 'SNOWMAN':
                    t_sh = self.topol_shorter_loop(bridge2,beg=bridge1[0],end=bridge1[1])
                elif topol_type == 'THETA':
                    t_sh = self.topol_theta_loop(bridge1,bridge2)
                elif topol_type == 'DEGENERATION' and bridge1[1]==bridge2[0]:
                    t_sh = self.topol_theta_loop(bridge1,bridge2)
                else:
                    continue
                t_in1 = self.topol_internal_loop(bridge1)
                t_in2 = self.topol_internal_loop(bridge2)
                with open('topology.txt', 'a+') as f:
                    f.write('  Shorter loop {:d}-{:d}_{:d}-{:d}: direct {}\n'.format(bridge1[0], bridge1[1], bridge2[0], bridge2[1], t_sh))
                    f.write('  Internal loop 1 {:d}-{:d}: direct {}\n'.format(bridge1[0], bridge1[1], t_in1))
                    f.write('  Internal loop 2 {:d}-{:d}: direct {}\n'.format(bridge2[0], bridge2[1], t_in2))
                    f.write('=====================\n')
                constituent_knots['in1_{:d}-{:d}'.format(bridge1[0],bridge1[1])] = t_in1
                constituent_knots['in2_{:d}-{:d}'.format(bridge2[0],bridge2[1])] = t_in2
                constituent_knots['sh_{:d}-{:d}_{:d}-{:d}'.format(bridge1[0],bridge1[1],bridge2[0],bridge2[1])] = t_sh
        return constituent_knots
        
    def topol_default_loop(self):
        coords = [self.coords[i][1:] for i in range(1,self.end_ndx+1)]

        with open('loop_xyzs/{}_default.xyz'.format(self.name), 'w') as f:
            f.write(self.format_coords(coords))
        return alexander(coords, max_cross=60, closure=2), alexander(coords, max_cross=60, closure=0), alexander(coords, max_cross=60, closure=1)

    def topol_internal_loop(self, bridge):
        coords = [self.coords[i][1:] for i in range(bridge[0],bridge[1]+1)]
        with open('loop_xyzs/{}_in_{:d}-{:d}.xyz'.format(self.name,bridge[0],bridge[1]), 'w') as f:
            f.write(self.format_coords(coords))
        ndx1, x1, y1, z1 = self.coords[bridge[0]]
        ndx2, x2, y2, z2 = self.coords[bridge[1]+1]
        dist = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        return alexander(coords, max_cross=60, closure=0)#, dist

    def topol_shorter_loop(self, bridge, beg=1, end=None):
        if end==None:
            end = self.end_ndx
        coords = [self.coords[i][1:] for i in range(beg,bridge[0]+1)]
        coords += [self.coords[i][1:] for i in range(bridge[1],end+1)]
        with open('loop_xyzs/{}_sh_{:d}-{:d}_b{:d}_e{:d}.xyz'.format(self.name,bridge[0],bridge[1],beg,end), 'w') as f:
            f.write(self.format_coords(coords))
        if beg == 1 and end == self.end_ndx:
            return alexander(coords, max_cross=60, closure=2), alexander(coords, max_cross=60, closure=0), alexander(coords, max_cross=60, closure=1)
        return alexander(coords, max_cross=60, closure=0)

    def topol_theta_loop(self, bridge1, bridge2):
        coords = [self.coords[i][1:] for i in range(bridge1[0],bridge2[0]+1)]
        coords += [self.coords[i][1:] for i in range(bridge2[1],bridge1[1]-1,-1)]
        with open('loop_xyzs/{}_th_{:d}-{:d}-{:d}-{:d}.xyz'.format(self.name,bridge1[0],bridge2[0],bridge2[1],bridge1[1]), 'w') as f:
            f.write(self.format_coords(coords))
        return alexander(coords, max_cross=60, closure=0)


if __name__ == '__main__':
    #auk = Protein('1aukA_CA.xyz', [12,263,264], [(282,396)])
    #ei6_1 = Protein('1ei6A_CA.xyz', [24,63,239,240,241],[(82,325),(22,57)])
    #ei6_2 = Protein('1ei6A_CA.xyz', [201,205,367],[(82,325),(22,57)])
    ejj_1 = Protein('1ejjA_CA.xyz', [10,60,442,443])
    ejj_2 = Protein('1ejjA_CA.xyz', [401,405,460])
    #protein1 = Protein('2k0aA_CA.xyz', [14, 49, 52, 89])
    #protein2 = Protein('2k0aA_CA.xyz', [26, 29, 61, 64])
    #protein3 = Protein('2k0aA_CA.xyz', [33, 36, 76, 79])