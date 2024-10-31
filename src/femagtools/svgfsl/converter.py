import sys
import pathlib
import argparse
import logging
import logging.config
from femagtools.dxfsl.converter import convert

logger = logging.getLogger(__name__)

def main():
    argparser = argparse.ArgumentParser(
        description='Process SVG file and create a plot or FSL file.')
    argparser.add_argument('svgfile',
                           help='name of SVG file')
    argparser.add_argument('--rotor',
                           help='rotor without airgap in/out',
                           dest='rotor',
                           default='')
    argparser.add_argument('--stator',
                           help='stator without airgap in/out',
                           dest='stator',
                           default='')
    argparser.add_argument('--EESM',
                           help='Electric Excited Synchronous Motor',
                           dest='EESM',
                           action="store_true",
                           default=False)
    argparser.add_argument('-p', '--plot',
                           help='show plots',
                           dest='show_plots',
                           action="store_true")
    argparser.add_argument('--areas',
                           help='show all areas',
                           dest='show_areas',
                           action="store_true")
    argparser.add_argument('--id',
                           help='write id of areas',
                           dest='write_id',
                           action="store_true")
    argparser.add_argument('-f', '--fsl',
                           help='create fsl',
                           dest='write_fsl',
                           action="store_true")
    argparser.add_argument('-v', '--view',
                           help='show a view only',
                           dest='view',
                           action="store_true")
    args = argparser.parse_args()

    loglevel = logging.INFO
    logging.basicConfig(level=loglevel,
                        format='%(asctime)s %(message)s',
                        filemode='w')
    part = ()
    if args.stator:
        if args.rotor:
            logger.error("Stator or Rotor expected")
            sys.exit(1)
        part = ('stator', args.stator)
    elif args.rotor:
        part = ('rotor', args.rotor)

    if not args.write_fsl:
        if not (args.show_plots or args.show_areas or args.view):
            args.write_fsl = True

    res = convert(args.svgfile,  # SVG-Filename
                  part=part,
                  EESM=args.EESM,
                  view_only=args.view,
                  show_plots=args.show_plots,
                  show_areas=args.show_areas,
                  write_fsl=args.write_fsl,
                  write_id=args.write_id)

    if args.write_fsl:
        if res is not None:
            p = pathlib.Path(args.svgfile)
            basename = p.stem
            pathlib.Path(basename + '.fsl').write_text('\n'.join(res['fsl']))

if __name__ == '__main__':
    main()
