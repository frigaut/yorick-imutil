%define name yorick-imutil
%define version 0.5.7
%define release gemini2008apr29

Summary: yorick library for image manipulation
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{version}.tar.bz2
License: BSD
Packager: Francois Rigaut <frigaut@gemini.edu>
Group: Development/Languages
Url: http://www.maumae.net/yorick/doc/plugins.php
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Requires: yorick >= 2.1 yorick-yutils


%description
Compiled routines for basic but fast image manipulation. Includes 2d
bilinear and spline2 interpolation, clipping, 2d dist generator,
binning, image rotation, cartesian to polar coordinate transform,
gaussian and poisson random generator, fast sort and fast median. All
of these functions, with the exceptions of spline2, exist in yorick or
the yutils package, but these versions are 2 to 10x faster, being
specialized for 2d arrays (hence the name imutil). This plugin is
64bits safe.

%prep
%setup -q

%build
yorick -batch make.i
make
if [ -f check.i ] ; then
   mv check.i %{name}_check.i
fi;
  
%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/lib
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/i0
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/i
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/i-start
mkdir -p $RPM_BUILD_ROOT/usr/share/doc/yorick-imutil
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/packages/installed

install -m 755 imutil.so $RPM_BUILD_ROOT/usr/lib/yorick/lib
install -m 644 *.i $RPM_BUILD_ROOT/usr/lib/yorick/i0
install -m 644 %{name}_check.i $RPM_BUILD_ROOT/usr/lib/yorick/i
install -m 644 *_start.i $RPM_BUILD_ROOT/usr/lib/yorick/i-start
install -m 644 LICENSE $RPM_BUILD_ROOT/usr/share/doc/yorick-imutil
install -m 644 README $RPM_BUILD_ROOT/usr/share/doc/yorick-imutil
install -m 644 imutil.info $RPM_BUILD_ROOT/usr/lib/yorick/packages/installed

rm $RPM_BUILD_ROOT/usr/lib/yorick/i0/*check.i
rm $RPM_BUILD_ROOT/usr/lib/yorick/i0/*_start.i


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
/usr/lib/yorick/lib/imutil.so
/usr/lib/yorick/i0/*.i
/usr/lib/yorick/i/%{name}_check.i
/usr/lib/yorick/i-start/*_start.i
/usr/share/doc/yorick-imutil/*
/usr/lib/yorick/packages/installed/*

%changelog
* Tue Jan 09 2008 <frigaut@users.sourceforge.net>
- included the info file for compat with pkg_mngr

* Mon Dec 31 2007 <frigaut@users.sourceforge.net>
- new distro directory structure
- updated cvs
