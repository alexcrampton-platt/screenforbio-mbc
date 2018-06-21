#!/usr/bin/perl

$maxlevel=shift;

while (<>) {
    print if ((split)[2] <= $maxlevel);
}
