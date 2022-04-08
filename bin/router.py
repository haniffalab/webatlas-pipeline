#!/usr/bin/env python3

import fire

def main(**args):
    type = args.pop('type')
    print(type)
    print(args)

if __name__ == "__main__":
    fire.Fire(main)