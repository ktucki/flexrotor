#!/usr/bin/env python3
"""Convert a myRIO/LabVIEW soil bin record (.xlsx, pairs of timestamp/value columns)
to the reduced CSV read by FlexRotor:
    t_us, shaft_cnt, disc_cnt, cart_cnt, torque_x10, speed_x10
Column letters of the source sheet can be changed with options.
Usage: python myrio_xlsx2csv.py record.xlsx [out.csv]
Requires: openpyxl
"""
import argparse, csv, os, sys

def col(letter):
    n = 0
    for ch in letter.upper():
        n = n * 26 + ord(ch) - 64
    return n - 1

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('xlsx'); ap.add_argument('csv', nargs='?')
    ap.add_argument('--cols', default='A,B,D,F,H,J', help='source columns: time,shaft,disc,cart,torque,speed')
    a = ap.parse_args()
    from openpyxl import load_workbook
    idx = [col(c) for c in a.cols.split(',')]
    out = a.csv or os.path.splitext(a.xlsx)[0] + '.csv'
    ws = load_workbook(a.xlsx, read_only=True, data_only=True).worksheets[0]
    n = 0
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t_us', 'shaft_cnt', 'disc_cnt', 'cart_cnt', 'torque_x10', 'speed_x10'])
        for row in ws.iter_rows(values_only=True):
            try:
                vals = [row[i] for i in idx]
            except IndexError:
                continue
            if not all(isinstance(v, (int, float)) for v in vals):
                continue
            w.writerow([int(v) if float(v).is_integer() else v for v in vals]); n += 1
    print(f'{out}: {n} samples', file=sys.stderr)

if __name__ == '__main__':
    main()
