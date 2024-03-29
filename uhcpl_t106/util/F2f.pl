#!/usr/bin/perl

@dirs = qw (modules src);

foreach $dir (@dirs) {

    chdir $dir;

    opendir DIR, ".";
    @files = grep /\.f90$/, readdir DIR;
    closedir DIR;

    foreach $file (@files) {
	$new = $file;
	$new =~ s/\.f90/\.F90/;
	rename $file, $new;
    }

    change_makefile ();

    chdir "..";
}

change_configure ();

exit;

sub change_configure {

    my $conffile = "configure";
    my $workfile = $conffile."tmp";

    open CONF, "<$conffile" or die "Cannot open configure!\n"; 
    open TMP,  ">$workfile" or die "Cannot open tmpfile!\n";

    while (<CONF>) {
	s/master\.f90/master\.F90/;
	s/for ac_prog in frt f90 xlf90/for ac_prog in f95 frt f90 xlf90/;
	print TMP $_;
    }

    close CONF;
    close TMP;

    rename $conffile, "configure.bak";
    rename $workfile, $conffile;
}

sub change_makefile {

    my $makefile = "Makefile";
    my $workfile = $makefile."tmp";

    open MAKE, "<$makefile" or die "Cannot open Makefile!\n"; 
    open TMP,  ">$workfile" or die "Cannot open tmpfile!\n";

    while (<MAKE>) {
	s/\.f90/\.F90/g;
	print TMP $_;
    }

    close MAKE;
    close TMP;

    rename $makefile, "Makefile.bak";
    rename $workfile, $makefile;
}
