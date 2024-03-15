#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import argparse

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-x", "--expect")
    return parser

def plot_coverage(df, baseline = None):
    plt.stem(df['coverage'])
    plt.xticks(range(1, len(df.index)+1))
    if baseline is not None:
        plt.axhline(y = float(baseline), color = 'r', linestyle='dashed')
    return plt


def main():
    args = build_parser().parse_args()
    df = pd.read_table(args.input).sort_values(by='sample')

    plot = plot_coverage(df, baseline = args.expect)
    plt.savefig(args.output)

if __name__ == '__main__':
    main()