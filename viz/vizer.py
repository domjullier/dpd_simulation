#!/usr/bin/env python2
# -*- coding: utf8 -*-

# Argumenty: (??)
# repeat - powtarzanie po dobiegnieciu do konca


# Program arguments.
import argparse, sys

aparser = argparse.ArgumentParser()
aparser.add_argument('result_dir',\
        help='directory containing visualisation data')
aparser.add_argument('--first', default=0, type=int,\
        help='first step of an animation (default: 0)')
aparser.add_argument('--last', default=sys.maxint, type=int,\
        help='last step of an animation (default: last)')
aparser.add_argument('--info_every', type=int, default=sys.maxint,\
        help='prints information every X steps.')
        

args = aparser.parse_args()

import os, logging
from wx import Timer
from mayavi import mlab

# Logging configuration.
log = logging.getLogger('vizer')
log.setLevel(logging.INFO)

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
    curStep = 0
    endStep = 0
    headFile = ''
    fileList = []

    def __init__(self, resDir, stepA, stepZ, infoEvery):
        log.debug('ResultProcessor starting.')
        self.resDir = resDir
        self.fileList = [os.path.join(resDir, f) for f in os.listdir(resDir)]
        self.fileList.sort()
        self.headFile = self.fileList.pop(0)

        self.curStep = stepA
        self.endStep = (len(self.fileList)-1, stepZ)[stepZ < len(self.fileList)]

        log.info('Results directory: %s, steps from %d to %d.',\
                self.resDir, self.curStep, self.endStep)

        self.infoEvery = infoEvery
        self.infoOffst = self.curStep % self.infoEvery

        log.debug('Animating files from %s, starting with %s up to %s.',\
                self.resDir, self.fileList[self.curStep],\
                self.fileList[self.endStep])

    def processHeader(self):
        log.info('Reading header from file %s.', self.headFile)
        with file(self.headFile, 'r') as infile:
            lines = [l.strip() for l in infile.xreadlines()]
            # TODO komentarze?
            #lines = [line for line in [l.strip() for l in infile.xreadlines()]\
            #        if line and not line.startswith('#')]

            # 1st line - world description
            desc = lines[0].split(',')

            log.debug('Headers 1st line: %s', lines[0])
            assert len(desc) > 4, 'Headers 1st line does not contain necessary\
                    fields!'

            self.wDesc.particles_n      = int(desc[1])
            self.wDesc.states_n         = int(desc[2])
            self.wDesc.states_interval  = float(desc[3])
            self.wDesc.world_x_size     = int(desc[4])
            self.wDesc.world_y_size     = int(desc[4])
            self.wDesc.world_z_size     = int(desc[4])

            # 2nd line empty, others - particle types
            for i in xrange(2,len(lines)):
                typel = lines[i].split(',')

                log.debug('Particle type: %s', lines[i])
                assert len(typel) > 2, 'Particle type description does not\
                        contain necessary fields!'

                parttype = ParticleType()
                parttype.n = int(typel[0])
                parttype.mass = float(typel[1])
                parttype.radius = float(typel[2])

                self.wDesc.ptypes.append(parttype)

    def processStep(self):
        if self.curStep > self.endStep:
            raise StopIteration

        x, y, z, s = [], [], [], []

        stepFile = self.fileList[self.curStep]

        if self.curStep % self.infoEvery == self.infoOffst:
            log.info('Animated step %d (file: %s).', self.curStep, stepFile)

        self.curStep += 1

        with file(stepFile, 'r') as infile:
            for line in infile.xreadlines():
                spl = line.strip().split(',')
                x.append(float(spl[1]))
                y.append(float(spl[2]))
                z.append(float(spl[3]))
                s.append(float(spl[0]))

            assert len(x) == self.wDesc.particles_n
            assert len(y) == self.wDesc.particles_n
            assert len(z) == self.wDesc.particles_n
            assert len(s) == self.wDesc.particles_n
            return (x,y,z,s)


class Animator(Timer):

    def __init__(self, mlab_source, fileproc):
        Timer.__init__(self)
        self.ms = mlab_source
        self.fproc = fileproc
        self.stepDelay = 1; #ms

    def Animate(self):
        log.debug('Animator starting.')
        self.Start(self.stepDelay, True)

    def Notify(self):
        try:
            x,y,z,s = self.fproc.processStep()
            self.ms.set(x = x, y = y, z = z, scalar=s)
            self.Start(self.stepDelay, True)
        except StopIteration:
            log.info('Reached the end of animation.')

def main():
    log.info('Program starting.')

    rp = ResultProcessor(args.result_dir, args.first, args.last,\
            args.info_every)

    rp.processHeader()
    s = rp.processStep()

    s = mlab.points3d(s[0], s[1], s[2], s[3], colormap="autumn",
            scale_factor = rp.wDesc.world_x_size / 100.)
            #scale_factor=10)

    anim = Animator(s.mlab_source, rp)
    anim.Animate()

    mlab.show()

if __name__ == '__main__':
    main()

