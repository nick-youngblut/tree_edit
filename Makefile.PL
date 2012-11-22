use 5.006;
use strict;
use warnings FATAL => 'all';
use inc::Module::Install;

name     'tree_edit';
all_from 'lib/tree_edit.pm';
author   q{Nick Youngblut <nyoungb2@illinois.edu>};
license  'perl';

perl_version 5.006;

tests_recursive('t');

resources (
   #homepage   => 'http://yourwebsitehere.com',
   #IRC        => 'irc://irc.perl.org/#tree_edit',
   license    => 'http://dev.perl.org/licenses/',
   #repository => 'git://github.com/Nick Youngblut/tree_edit.git',
   bugtracker => 'http://rt.cpan.org/NoAuth/Bugs.html?Dist=tree_edit',
);

configure_requires (
   'Module::Install' => 0,
);

build_requires (
   'Test::More' => 0,
);

requires (
   #'ABC'              => 1.6,
   #'Foo::Bar::Module' => 5.0401,
);

install_as_cpan;
auto_install;
WriteAll;
