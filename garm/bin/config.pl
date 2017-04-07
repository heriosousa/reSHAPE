#!/usr/bin/env perl
### CHECK DEPENDENCIES FOR GARM ###

use strict;
use warnings;
use Getopt::Std;
use File::Temp;

my $shell = $ENV{SHELL};
my $sh = (split /\//, $shell)[$#_];
print "Your shell is $shell\nConfiguring GARM env variables\n"; #<STDIN>;
my $pwd = $ENV{PWD};
my $garmbin = $ENV{GARMBIN};
$garmbin = "$pwd/bin" unless $garmbin;
my $mumbin = $ENV{MUMBIN};
my $amosbin = $ENV{AMOSBIN};
my $garmlib = $ENV{GARMLIB};
$garmlib = "$pwd/lib" unless $garmlib;

my @amosapps = qw(toAmos bank-transact nucmer2ovl sort2 tigger make-consensus bank2contig bank2fasta listReadPlacedStatus dumpreads);
my $totalerror = 0;
my %amospaths;
my $already = 0;
if ($garmbin and $garmlib and $mumbin and $amosbin) {
	$already = 1;
	print "### Your enviromental variables are already defined:\n";
	print "$garmbin\n$garmlib\n$mumbin\n$amosbin\n\n";
	print "### Checking apps\n";
	my $nucmer = "$mumbin/nucmer";
	$totalerror = CHECKAP("$nucmer");
	foreach my $app (@amosapps) {
		$totalerror += CHECKAP("$amosbin/$app");
	}
} elsif ($garmbin and $garmlib)  {
	print "GARM enviromental variables defined:\n";
	if ($mumbin) {
		print "MUMmer path defined:\n$mumbin\n";
	} else {
		print "MUMmer path NOT defined! Guessing path from applications:\n";
		my $nucmer = `which nucmer`;
		chomp $nucmer;
		$totalerror = CHECKAP("$nucmer");
		my @path = split(/\//, $nucmer);
		pop @path;
		my $newpath = join("/", @path);
		$mumbin = $newpath;
	}
	if ($amosbin) {
		print "AMOS path defined:\n$amosbin\n";
		print "Checking apps\n"; <STDIN>;
		foreach my $app (@amosapps) {
			$totalerror += CHECKAP("$amosbin/$app");
		}
	} else {
		print "AMOS path NOT defined! Guessing path from applications:\n";
		my $toAmos = `which toAmos`;
		my $error = `echo $?`;
		my $banktran = `which bank-transact`;
		$error += `echo $?`;
		my $nucovl = `which nucmer2ovl`;
		$error += `echo $?`;
		my $sort2 = `which sort2`;
		$error += `echo $?`;
		my $ovlOVL = `which ovl2OVL`;
		$error += `echo $?`;
		my $tigger = `which tigger`;
		$error += `echo $?`;
		my $makecon = `which make-consensus`;
		$error += `echo $?`;
		my $bank2con = `which bank2contig`;
		$error += `echo $?`;
		my $bank2fas = `which bank2fasta`;
		$error += `echo $?`;
		my $listreads = `which listReadPlacedStatus`;
		$error += `echo $?`;
		my $dumpreads = `which dumpreads`;
		chomp $error;
		if ($error) {
			print "$error application(s) were not found in your path\n";
			print "Please install AMOS package and be sure that is in your path or define it!\n";
			$totalerror = $error;
		}
		my @guess;
		push (@guess,$toAmos );
		push (@guess,$banktran);
		push (@guess,$nucovl);
		push (@guess,$sort2);
		push (@guess,$ovlOVL);
		push (@guess,$tigger);
		push (@guess,$makecon);
		push (@guess,$bank2con);
		push (@guess,$bank2fas);
		push (@guess,$listreads);
		push (@guess,$dumpreads);
		print "Checking apps\n";
		foreach my $app (@guess) {
			chomp $app;
			my @path = split(/\//, $app);
			pop @path;
			my $newpath = join("/", @path);
			#print $newpath; <STDIN>;
			$amospaths{$newpath}++;
			$totalerror += CHECKAP($app);
		}
	}
}

if ($totalerror) {
	die "GARM cannot be run due to several dependencies are missing\nPlease check the README.txt and install everything properly\n";
} else {
	my $checkpath = scalar (keys %amospaths);
	if ($checkpath and $checkpath == 1) {
		$amosbin = (keys %amospaths)[0]; 
	} elsif ( $checkpath > 1) {
		die "The AMOS applications are distributed in several locations...\nPlease, install or link those applications in just one place 
and make it available in your path\n";
	}
	open SH, ">$pwd/GARMenv.sh" || die "Cannot open GARMenv.sh: $!\n";
	if ($shell =~ /tcsh/) {
		print SH "setenv GARMBIN $pwd/bin\n"; #<STDIN>;
		print SH "setenv GARMLIB $pwd/lib\n"; #<STDIN>;
		print SH "setenv MUMBIN $mumbin\n";
		print SH "setenv AMOSBIN $amosbin\n";
	}
	elsif ($shell =~ /bash/) {
		print SH "export GARMBIN=$pwd/bin\n";
		print SH "export GARMLIB=$pwd/lib\n";
		print SH "export MUMBIN=$mumbin\n";
		print SH "export AMOS=$amosbin\n";
	}
	print "\nYou can run now GARM without problems :D (hopefully...)\n";
	unless ($already) {
		print "You can modify to your .cshrc or .bashrc file and add the variables in the file GARMenv.sh\n";
		print "or you have to source it BEFORE running GARM every time\n(source GARMenv.sh)\n";
	}
	close SH;
}

sub CHECKAP {
	my $app = shift;
    my $help;
	print "Checking $app\n"; #<STDIN>;
    if( $shell =~ /tcsh/ ){
        $help = `$app -h >&/dev/null`;
    }else{
        $help = `$app -h 2>&1`;
    }
	my $error = `echo $?`;
	chomp $error;
	#print $version;
	if ($error) {
		if ($help =~ /Command not found/) {
			print "GARM can't run if $app is not installed or running properly\nPlease install and check it first!\n";
		} else {
			print "OK\n";
			$error = 0;
		}
	} else {
		print "OK\n";
	}

return ($error);

}
