"""
  femagtools.dxfsl.conv

  entry point for DXF to FSL conversion

"""
import sys
import os
import io
import femagtools
from femagtools.dxfsl.converter import convert
import argparse
import logging
import logging.config

logger = logging.getLogger(__name__)


def main():
    argparser = argparse.ArgumentParser(
        description='Process DXF file and create a plot or FSL file.')
    argparser.add_argument('dxfile',
                           help='name of DXF file')
    argparser.add_argument('--inner',
                           help='name of inner element',
                           dest='inner',
                           default='inner')
    argparser.add_argument('--outer',
                           help='name of outer element',
                           dest='outer',
                           default='outer')
    argparser.add_argument('--rotor',
                           help='rotor without airgap in/out',
                           dest='rotor',
                           default='')
    argparser.add_argument('--stator',
                           help='stator without airgap in/out',
                           dest='stator',
                           default='')
    argparser.add_argument('-a', '--airgap',
                           help='correct airgap',
                           dest='airgap',
                           type=float,
                           default=0.0)
    argparser.add_argument('--airgap2',
                           help='correct airgap',
                           dest='airgap2',
                           type=float,
                           default=0.0)
    argparser.add_argument('-t', '--symtol',
                           help='absolut tolerance to find symmetry axis',
                           dest='sym_tolerance',
                           type=float,
                           default=0.001)
    argparser.add_argument('--mindist',
                           help='minimal distance of spline control-points',
                           dest='mindist',
                           type=float,
                           default=0.01)
    argparser.add_argument('--rtol',
                           help='relative tolerance (pickdist)',
                           dest='rtol',
                           type=float,
                           default=1e-03)
    argparser.add_argument('--atol',
                           help='absolut tolerance (pickdist)',
                           dest='atol',
                           type=float,
                           default=0.005)
    argparser.add_argument('--da',
                           help='distance airgap',
                           dest='da',
                           type=float,
                           default=0.0)
    argparser.add_argument('--dy',
                           help='distance yoke',
                           dest='dy',
                           type=float,
                           default=0.0)
    argparser.add_argument('-s', '--split',
                           help='split intersections',
                           dest='split',
                           action="store_true")
    argparser.add_argument('-p', '--plot',
                           help='show plots',
                           dest='show_plots',
                           action="store_true")
    argparser.add_argument('--areas',
                           help='show all areas',
                           dest='show_areas',
                           action="store_true")
    argparser.add_argument('-f', '--fsl',
                           help='create fsl',
                           dest='write_fsl',
                           action="store_true")
    argparser.add_argument('-v', '--view',
                           help='show a view only',
                           dest='view',
                           action="store_true")
    argparser.add_argument('-k', '--korr',
                           help='show a view with korrections',
                           dest='view_korr',
                           action="store_true")
    argparser.add_argument('--png',
                           help='write png-file only',
                           dest='write_png',
                           action="store_true")
    argparser.add_argument('-d', '--debug',
                           help='print debug information in logfile',
                           dest='debug',
                           action="store_true")
    argparser.add_argument('-l', '--log',
                           help='print information in logfile',
                           dest='debug',
                           action="store_true")
    argparser.add_argument('--version',
                           help='show version of some packages',
                           dest='version',
                           action="store_true")
    argparser.add_argument('--debugger',
                           help='print debug information in logfile',
                           dest='debugger',
                           action="store_true")

    args = argparser.parse_args()

    logfilename = None
    loglevel = logging.INFO
    if args.debug:
        loglevel = logging.DEBUG
        logfilename = 'debugger.log'
        print("see log-messages in {}".format(logfilename))

    logging.basicConfig(level=loglevel,
                        format='%(asctime)s %(message)s',
                        filename=logfilename,
                        filemode='w')

    if args.version:
        logger.info("femagtools version: %s", femagtools.__version__)
        try:
            import networkx as nx
            logger.info("networkx version: %s", nx.__version__)
        except ImportError:  # ModuleNotFoundError:
            logger.info("networkx version: <networkx not available>")
        try:
            import matplotlib
            logger.info("matplotlib version: %s", matplotlib.__version__)
        except ImportError:  # ModuleNotFoundError:
            logger.info("matplotlib version: <matplotlib not available>")
        logger.info("Python: %s", sys.version)
        sys.exit(0)

    if args.airgap > 0.0:
        if args.airgap2 > 0.0:
            logger.info("Airgap is set from {} to {}"
                        .format(args.airgap, args.airgap2))
        else:
            logger.info("Airgap is set to {}".format(args.airgap))

    part = ()
    if args.stator:
        if args.rotor:
            logger.error("Stator or Rotor expected")
            sys.exit(1)
        part = ('stator', args.stator)
    elif args.rotor:
        part = ('rotor', args.rotor)
    if part:
        if args.airgap:
            logger.info('airgap in stator or rotor not possible')
            sys.exit(1)
        args.airgap = -1  # no airgap
        if part[1] not in ('in', 'out'):
            logger.info('{} has to be defined in/out'.format(part[0]))
            sys.exit(1)

    if not args.write_fsl:
        if not (args.show_plots or args.show_areas or args.view):
            args.write_fsl = True

    res = convert(args.dxfile,  # DXF-Filename
                  rtol=args.rtol,    # relative pickdist toleranz
                  atol=args.atol,    # absolute pickdist toleranz
                  symtol=args.sym_tolerance,
                  mindist=args.mindist,
                  split=args.split,
                  inner_name=args.inner,
                  outer_name=args.outer,
                  part=part,
                  airgap=args.airgap,
                  airgap2=args.airgap2,
                  da=args.da,  # distance airgap
                  dy=args.dy,  # distance yoke
                  view_only=args.view,
                  view_korr=args.view_korr,
                  show_plots=args.show_plots,
                  show_areas=args.show_areas,
                  write_fsl=args.write_fsl,
                  write_png=args.write_png,
                  debug_mode=args.debugger)

    if args.write_fsl:
        if res is not None:
            basename = os.path.basename(args.dxfile).split('.')[0]
            with io.open(basename + '.fsl', 'w', encoding='utf-8') as f:
                f.write('\n'.join(res['fsl']))


if __name__ == "__main__":
    loglevel = logging.INFO

    main()
