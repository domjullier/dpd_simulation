#!/usr/bin/env python2
# -*- coding: utf8 -*-

# Argumenty: (??)
# repeat - powtarzanie po dobiegnieciu do konca
# first, last - kroki od-do.
# logowanie jakiejś informacji co x kroków?

import argparse, os, logging
from wx import Timer
from mayavi import mlab

# Program arguments.
aparser = argparse.ArgumentParser()
aparser.add_argument('result_dir', help='directory containing visualisation\
        data')

args = aparser.parse_args()

# Logging configuration.
log = logging.getLogger('vizer')
log.setLevel(logging.DEBUG)

han = logging.StreamHandler()
han.setLevel(logging.DEBUG)

log.addHandler(han)



class ParticleType(object):
    n = 0
    mass = 0
    radius = 0

class WorldDescription(object):
    particles_n = 0
    states_n = 0
    states_interval = 0
    world_x_size = 0
    world_y_size = 0
    world_z_size = 0
    ptypes = []

class ResultProcessor(object):
    resDir = ''
    wDesc = WorldDescription()

    def __init__(self, resDir):
        self.resDir = resDir

    def processHeader(self):
        with file(os.path.join(self.resDir, 'result'), 'r') as infile:
            lines = [l.strip() for l in infile.xreadlines()]
            #lines = [l.strip() for l in infile.xreadlines() if l and l.strip()[0] != '#']

            # 1st line - world description
            desc = lines[0].split(',')

            self.wDesc.particles_n      = int(desc[1])
            self.wDesc.states_n         = int(desc[2])
            self.wDesc.states_interval  = float(desc[3])
            self.wDesc.world_x_size     = int(desc[4])
            self.wDesc.world_y_size     = int(desc[4])
            self.wDesc.world_z_size     = int(desc[4])

            # 2nd line empty, others - particle types
            for i in xrange(2,len(lines)):
                typel = lines[i].split(',')
                log.debug('Typeline[%d]: %s', i, typel[0])

                parttype = ParticleType()
                parttype.n = int(typel[0])
                parttype.mass = float(typel[1])
                parttype.radius = float(typel[2])

                self.wDesc.ptypes.append(parttype)

    def processStep(self, n):
        x, y, z, s = [], [], [], []
        with file(os.path.join(self.resDir, 'result_step_' + str(n)), 'r') as infile:
            for line in infile.xreadlines():
                spl = line.strip().split(',')
                x.append(float(spl[1]))
                y.append(float(spl[2]))
                z.append(float(spl[3]))
                s.append(float(spl[0]))

            assert len(x) == self.wDesc.particles_n # TODO dodać kolejne??
            return (x,y,z,s)


class Animator(Timer):

    def __init__(self, mlab_source, fileproc):
        Timer.__init__(self)
        self.ms = mlab_source
        self.fproc = fileproc
        self.step_n = 1
        self.step_delay = 1; #ms

    def Animate(self):
        log.debug('Animator starting.')
        self.Start(self.step_delay, True)

    def Notify(self):
        try:
            s = self.fproc.processStep(self.step_n)
            self.step_n += 1
            self.ms.set(x = s[0], y = s[1], z = s[2], scalar=s[3])
            self.Start(self.step_delay, True)
        except EOFError:
            log.warn('End of file.') # TODO na pewno?
        except IOError:
            log.info('Reached end of animation.')

def main():
    log.info('Program starting.')
    log.info('Results directory: %s', args.result_dir)
    log.info('Steps: %d to %d.', 0, len(os.listdir(args.result_dir))-2)

    rp = ResultProcessor(args.result_dir)
    rp.processHeader()
    s = rp.processStep(0)

    s = mlab.points3d(s[0], s[1], s[2], s[3], colormap="autumn",
            scale_factor=15)

    anim = Animator(s.mlab_source, rp)
    anim.Animate()

    mlab.show()

if __name__ == '__main__':
    main()

