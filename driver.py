import os, argparse, pickle

from actsims.simTools import freqsInPsas, getActpolNoiseSim, getActpolCmbFgSim
from actsims.flipperDict import flipperDict

from pixell import enmap

parser = argparse.ArgumentParser(description="Generate simulations for act")
parser.add_argument("--iteration", type=int, default=0, help="iteration number")
parser.add_argument("--set", type=int, default=0, help="cmb set number")
parser.add_argument("--patch", default='deep5', help="patch")
parser.add_argument("--season", default='s13', help="season")
parser.add_argument("--array", default='pa1', help="array")
parser.add_argument("--noise", help="noise parameter file, i.e. templateInputsMr3c", required=True)
parser.add_argument("--signal", help="signal parameter file, i.e. signal.dict", required=True)
parser.add_argument("--verbose", action="store_true", default=False, help="verbose option")
parser.add_argument("--sim_type", default="all", help="type of data to simulate")
parser.add_argument("--window", action="store_true", default=True, help="apply window function")
parser.add_argument("--beam", action="store_true", default=True, help="do beam convolution")
parser.add_argument("--diagonal", action="store_true", default=False, help="diagonal noise only")
parser.add_argument("--split", type=int, default=0, help="number of splits wanted")
parser.add_argument("--modulation", action="store_true", default=True, help="apply modulation")
parser.add_argument("--cmb_type", default="LensedCMB", help="type of the cmb data of interests")

args = parser.parse_args()

##############
# parameters #
##############

# output options
output_dir = "outputs"
prefix = "test"

# random seeds
cmbSeedInd = 0
fgSeedInd = 1
phiSeedInd = 2
noiseSeedInd = 3

################
# main program #
################

# load parameter files for noise simulation
nDict = flipperDict()
nDict.read_from_file(args.noise)

# load parameter files for signal simulation
sDict = flipperDict()
sDict.read_from_file(args.signal)

# patch season and array information
psa = '%s_%s_%s' % (args.patch, args.season, args.array)

# figure out what frequencies correspond to this array, using the
# function defined above
psaFreqs = freqsInPsas(psa, nDict['freqsInArrays'])

# unroll psa names (normally stored as a nested  list of lists)
psaList = [item for sublist in nDict['psaList'] for item in sublist]

# check if psa of interest is in the psa list
if psa not in psaList:
    raise ValueError('psa %s not found in psaList; options are ' % (psa ), psaList)

psa_id = psaList.index(psa)
noiseSeed = (args.set, psa_id, noiseSeedInd*4+args.split, args.iteration)

# load up one sample map, just to get the shape and wcs info.  Do this for "I" at one frequency
filename = os.path.join(nDict['dataMapDir'],
                        'totalWeightMapI_{psa}_{freq}_fromenlib.fits'.format(psa=psa,
                                                                             freq=psaFreqs[0]))
sampleMap = enmap.read_map(filename)


# Note: Foreground seed is the same for every sky patch, season, and
# frequency!  This is because they are used to generate fullsky alm's
foregroundSeed = (args.set, 0, fgSeedInd, args.iteration)


# check output folder
if not os.path.exists(output_dir):
    print("Warning: %s doesn't exist, creating one now" % output_dir)
    os.makedirs(output_dir)


# noise simulation
if args.sim_type == "all" or args.sim_type == "noise":
    results = getActpolNoiseSim(noiseSeed=noiseSeed, psa=psa,
                                noisePsdDir=nDict['dataMapDir'],
                                freqs=psaFreqs, verbose=args.verbose,
                                noiseDiagsOnly=args.diagonal,
                                splitWanted=args.split)

    filename = os.path.join(output_dir, "%s_noise.pickle" % prefix)
    with open(filename, "wb") as f:
        print("Saving file: %s..." % filename)
        pickle.dump(results, f)


# cmb simulation
if args.sim_type == 'all' or args.sim_type == 'cmb':
    results = getActpolCmbFgSim(beamfileDict=sDict['beamNames'],
                                shape=sampleMap.shape, wcs=sampleMap.wcs,
                                iterationNum=args.iteration,
                                cmbDir=sDict['cmbDir'], freqs=psaFreqs,
                                psa=psa, cmbSet=args.set,
                                doBeam=args.beam,
                                applyWindow=args.window,
                                verbose=args.verbose,
                                cmbMaptype=args.cmb_type,
                                foregroundSeed=foregroundSeed,
                                simType='cmb',
                                foregroundPowerFile=sDict['foregroundPowerFile'],
                                applyModulation=args.modulation)

    filename = os.path.join(output_dir, "%s_cmb.pickle" % prefix)
    with open(filename, "wb") as f:
        print("Saving file: %s..." % filename)
        pickle.dump(results, f)


if args.sim_type == 'all' or args.sim_type == 'foregrounds':
    results = getActpolCmbFgSim(beamfileDict=sDict['beamNames'],
                                shape=sampleMap.shape, wcs=sampleMap.wcs,
                                iterationNum=args.iteration,
                                cmbDir=sDict['cmbDir'], freqs=psaFreqs,
                                psa=psa, cmbSet=args.set,
                                doBeam=args.beam,
                                applyWindow=args.window,
                                verbose=args.verbose,
                                cmbMaptype=args.cmb_type,
                                foregroundSeed=foregroundSeed,
                                simType='foregrounds',
                                foregroundPowerFile=sDict['foregroundPowerFile'],
                                applyModulation=args.modulation)

    filename = os.path.join(output_dir, "%s_foreground.pickle" % prefix)
    with open(filename, "wb") as f:
        print("Saving file: %s..." % filename)
        pickle.dump(results, f)


print("Done!")
