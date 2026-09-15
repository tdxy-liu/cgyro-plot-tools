"""Read-only installation and optional case check; never opens a GUI window."""
import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', help='optional CGYRO case directory')
    parser.add_argument('--verify-all', action='store_true',
                        help='decode/check every committed block (can take time)')
    args = parser.parse_args(argv)
    if args.verify_all and not args.case:
        parser.error('--verify-all requires --case')
    print('Python:', sys.executable)
    try:
        import cgyro_comparison_bootstrap as bootstrap
        import pygacode
        import cgyro_fullt_reader as codec
        print('pygacode:', pygacode.__file__)
        print('FULLT reader (bundled):', codec.__file__)
        print('Runtime ABI:', bootstrap._runtime_abi)
        zstd = codec._zstd()
        sample = b'CGYRO decoder health check\0' * 16
        restored = zstd.ZstdDecompressor().decompress(
            zstd.ZstdCompressor(level=3).compress(sample))
        if restored != sample:
            raise RuntimeError('zstandard round-trip failed')
        print('zstandard:', zstd.__version__, zstd.__file__)
        if args.case:
            from cgyro_fullt_compressed_bridge import open_fullt
            reader = open_fullt(args.case)
            if reader is None:
                raise RuntimeError('No selected FULL_T_ASYM data found')
            with reader:
                print('Selected data:', reader.path)
                print('Layout:', reader.layout)
                print('Available records:', reader.nt)
                count = 0
                if args.verify_all:
                    for t in range(reader.nt):
                        for q in range(reader.layout.ny):
                            reader.block(t, q)
                            count += 1
                elif reader.nt:
                    reader.block(reader.nt - 1, min(1, reader.layout.ny - 1))
                    count = 1
                print('Decoded blocks:', count,
                      '(FTZ checksums checked; native raw has no stored checksum)')
        print('OK: compressed-reading runtime is ready; no data files modified.')
        return 0
    except Exception as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
