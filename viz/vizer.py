#!/usr/bin/env python2
# -*- coding: utf8 -*-

# Argumenty: (??)
# repeat - powtarzanie po dobiegnieciu do konca

import argparse, random

aparser = argparse.ArgumentParser()
aparser.add_argument('infile', help='file containing particle data')

args = aparser.parse_args()

class WorldDescription(object):
    particles_n = 0
    states_n = 0
    states_interval = 0
    world_size = 0


class FileProcessor(object):

    def __init__(self, filename):
        self.ifile = open(filename, 'r')
        self.desc = WorldDescription()
        self.ptype = []

    def process_header(self):
        try:
            header = self.ifile.next().strip().split(',')

            self.desc.particles_n = int(header[1])
            self.desc.states_n = int(header[2])
            self.desc.states_interval = float(header[3])
            self.desc.world_size = int(header[4])

            assert(self.ifile.next().strip() == '')

            typedict = {}
            #for p in xrange(self.desc.particles_n):
            #    types = self.ifile.next().strip().split(',')
            #    if types[0] not in typedict:
            #        print 'Types[0]:', types[0]
            #        typedict[types[0]] = float(types[0])
            #    self.ptype.append(typedict[types[0]])
                    
            line = self.ifile.next().strip()
            while line != '':
                types = line.split(',')
                if types[0] not in typedict:
                    print 'Types[0]:', types[0]
                    typedict[types[0]] = float(types[0])
                self.ptype.append(typedict[types[0]])
                line = self.ifile.next().strip()

            #assert(self.ifile.next().strip() == '')
        except StopIteration:
            print 'Error - illformed input file.'
    
    def process_step(self):
        x, y, z, s = [], [], [], []
        try:
            for p in xrange(self.desc.particles_n):
                spl = self.ifile.next().strip().split(',')
                x.append(float(spl[1]))
                y.append(float(spl[2]))
                z.append(float(spl[3]))
                s.append(float(spl[0]))
            assert(self.ifile.next().strip() == '')
        except StopIteration:
            raise EOFError()

        return (x, y, z, s)




from wx import Timer

step_ms = 40

class Animator(Timer):

    def __init__(self, mlab_source, fileproc):
        Timer.__init__(self)
        self.ms = mlab_source
        self.fproc = fileproc

    def Notify(self):
        #print 'Notify!'
        try:
            s = self.fproc.process_step()
            #print s
            self.ms.set(x = s[0], y = s[1], z = s[2], scalar=s[3])
            self.Start(step_ms, True)
        except EOFError:
            print 'End of file.'



fp = FileProcessor(args.infile)
fp.process_header()

s = fp.process_step()

from mayavi import mlab
s = mlab.points3d(s[0], s[1], s[2], s[3], colormap="autumn", scale_factor=15)

anim = Animator(s.mlab_source, fp)
anim.Start(step_ms, True)



mlab.show()

