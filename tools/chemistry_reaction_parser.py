#!/usr/bin/env python3

import re
import pandas as pd
import argparse


p = argparse.ArgumentParser()
p.add_argument('csv_file', type=str, help='Input csv file')
p.add_argument('-convert_tex', type=str,
               help='Convert some LaTeX and save in new csv file')
p.add_argument('-default_length_unit', type=str, default='cm',
               help='Length unit to use for reactions')
p.add_argument('-comment', action='store_true',
               help='Include comments')
args = p.parse_args()

# re for a floating point number. The (?:  ) is a non-capturing group
number = r'[+-]?(?:[0-9]*[.])?[0-9]+'

# re for a float in scientific notation (e.g. 0.5e-10)
s_number = number + r'(?:[eEdD]' + number + ')?'


def reaction_to_regex(text):
    """Convert reaction text to a regular expression"""

    r = text

    # Check that coefficients are c1, c2, ..., cN
    c = ''.join(re.findall(r'[+-]?c([0-9])', r))
    if c not in '123456789':
        raise ValueError('Coefficients should appear in order c1, c2, ...')

    # Find +- signs in front of coefficients
    signs = re.findall(r'([+-]?)c[0-9]', r)
    c_signs = [int(s+'1') for s in signs]

    # Remove +- signs from coefficients
    r = re.sub(r'[+-](c[0-9])', r' \1', r)

    # Protect special symbols with []
    r = re.sub(r'\*', '[*]', r)
    r = re.sub(r'\+', '[+]', r)
    r = re.sub(r'\(', r'[(]', r)
    r = re.sub(r'\)', r'[)]', r)

    # Replace word boundaries and spaces by one or more spaces
    r = re.sub(r' ', r' *', r)
    r = re.sub(r'\b', r' *', r)

    # Convert c0, c1, ... to capture groups
    r = re.sub(r'c[0-9]', '(' + s_number + ')', r)

    # Allow for whitespace at beginning or end
    r = r'^\s*' + r + r'\s*$'

    return r, c_signs


formats = [
    r'c1',
    r'c1*(Td-c2)',
    r'c1*exp(-(c2/(c3+Td))**2)',
    r'c1*exp(-(Td/c2)**2)',
    r'c1*(300/Te)**c2',
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
]

rate_funcs = [
    {'text': x,
     'name': x.replace(' ', ''),
     'regex_signs': reaction_to_regex(x)
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
    [r'T_d', r'Td'],
    [r'T_e', r'Te'],
    [r'T_g', r'Tg'],
    # \to -> '->'
    [r'\\to', r'->'],
]


def replace_tex(text):
    for old, new in tex_replacements:
        text = re.sub(old, new, text)
    return text


reactions = pd.read_csv(args.csv_file, comment='#')

if args.convert_tex:
    reactions['reaction'] = [replace_tex(t) for t in reactions['reaction']]
    reactions['rate'] = [replace_tex(t) for t in reactions['rate']]
    reactions.to_csv(args.convert_tex, index=False)
else:
    for index, row in reactions.iterrows():
        r = row['reaction']
        f = row['rate']

        if 'length_unit' in row:
            length_unit = row['length_unit']
        else:
            length_unit = args.default_length_unit

        found_match = False

        for k in rate_funcs:
            m = re.match(k['regex_signs'][0], f.strip())
            if m:
                # Correct signs of coefficients
                signs = k['regex_signs'][1]
                coeffs = [float(x)*s for x, s in zip(m.groups(), signs)]
                coeffs = ' '.join([str(x) for x in coeffs])

                if args.comment and 'comment' in row:
                    print('# ' + row['comment'].strip())

                print(f'{r.strip()},{k["name"]},{coeffs},{length_unit}')
                found_match = True
                break
        if not found_match:
            print(f'** No match for {f.strip()}')
