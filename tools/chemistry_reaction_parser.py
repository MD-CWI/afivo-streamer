#!/usr/bin/env python3

import re
import pandas as pd
import argparse


p = argparse.ArgumentParser()
p.add_argument('csv_file', type=str, help='Input csv file')
p.add_argument('-convert_tex', action='store_true',
               help='Try to convert LaTeX reaction equations')
p.add_argument('-default_length_unit', type=str, default='cm',
               help='Length unit to use for reactions')
args = p.parse_args()

# re for a floating point number. The (?:  ) is a non-capturing group
number = r'[+-]?(?:[0-9]*[.])?[0-9]+'

# re for a float in scientific notation (e.g. 0.5e-10)
s_number = number + r'(?:[eEdD]' + number + ')?'


def reaction_to_regex(text):
    """Convert reaction text to a regular expression"""

    # Protect special symbols with []
    r = re.sub(r'\*', '[*]', text)
    r = re.sub(r'\+', '[+]', r)
    r = re.sub(r'\(', r'[(]', r)
    r = re.sub(r'\)', r'[)]', r)

    # Remove spaces in original description
    r = re.sub(r' ', r'', r)

    # Replace word boundaries by one or more spaces
    r = re.sub(r'\b', r' *', r)

    # Convert c0, c1, ... to capture groups
    r, count = re.subn(r'c[0-9]+', '(' + s_number + ')', r)

    # Allow for whitespace at beginning or end
    r = r'^\s*' + r + r'\s*$'
    return r


formats = [
    r'c1',
    r'c1*(Td-c2)',
    r'c1*exp(-(c2/(c3+Td))**2)',
    r'c1*exp(-(Td/c2)**2)',
    r'c1 * (300/Te)**c2',
    r'(c1*(kB_eV*Te+c2)**2-c3)*c4',
    r'c1*(Tg/300)**c2*exp(-c3/Tg)',
    r'c1*exp(-c2/Tg)',
    r'c1*Tg**c2',
    r'c1*(Tg/c2)**c3',
    r'c1*(300/Tg)**c2',
    r'c1*exp(-c2*Tg)',
    r'10**(c1+c2*(Tg-300))',
    r'c1*(300/Tg)**c2*exp(-c3/Tg)',
    r'c1*Tg**c2*exp(-c3/Tg)',
    r'c1*exp(-(c2/(c3+Td))**c4)',
    r'c1*exp(-(Td/c2)**c3)',
    r'c1*exp(-(c2/(kb*(Tg+Td/c3)))**c4)',
    # New reactions
    r'c1 * exp(c2 / Tg)',
]

rate_funcs = [
    {'text': x,
     'name': x.replace(' ', ''),
     'regex': reaction_to_regex(x)
     } for x in formats]

tex_replacements = [
    # 2.40\times10^{-7} -> 2.40e-7
    ['(' + number + r') *\\times *10\^\{(' + number + ')\}', r'\1e\2'],
    # x^{0.7} -> x**0.7
    [r'\^\{(' + number + r')\}', r'**\1'],
    # \frac{x/T} -> x/T
    [r'\\frac\{(' + number + ')\}\{(\w+)\}', r'\1/\2'],
    # \frac{T/x} -> T/x
    [r'\\frac\{(\w+)\}\{(' + number + ')\}', r'\1/\2'],
    # Add *
    [r'([0-9])\(', r'\1*('],
    # Add *
    [r'([0-9])\\', r'\1*\\'],
    # \exp -> exp
    [r'\\exp', r'exp'],
    # T_e -> Te
    [r'T_e', r'Te'],
    [r'T_g', r'Tg'],
    # \to -> '->'
    [r'\\to', r'->'],
]

reactions = pd.read_csv(args.csv_file, comment='#')

for index, row in reactions.iterrows():
    r = row['reaction']
    f = row['rate']

    if 'length_unit' in row:
        length_unit = row['length_unit']
    else:
        length_unit = args.default_length_unit

    if args.convert_tex:
        for old, new in tex_replacements:
            f = re.sub(old, new, f)
            r = re.sub(old, new, r)

    found_match = False

    for k in rate_funcs:
        m = re.match(k['regex'], f.strip())
        if m:
            coeffs = ' '.join(m.groups())

            if 'comment' in row:
                print('# ' + row['comment'].strip())

            print(f'{r.strip()},{k["name"]},{coeffs},{length_unit}')
            found_match = True
    if not found_match:
        print(f'** No match for {f.strip()}')
