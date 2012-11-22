#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'tree_edit' ) || print "Bail out!\n";
}

diag( "Testing tree_edit $tree_edit::VERSION, Perl $], $^X" );
