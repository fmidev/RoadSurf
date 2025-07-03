%define NAME roadsurf
Summary: Road surface weather model
Name: %{NAME}
Version: 25.7.3
Release: 1%{?dist}.fmi
License: MIT
Group: Development/Tools
URL: https://github.com/fmidev/RoadSurf
Source0: %{name}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot-%(%{__id_u} -roadn)
BuildRequires: rpm-build
BuildRequires: gcc-gfortran
BuildRequires: make
Provides: RoadSurf

%description
Road surface weather model developed by the Finnish Meteorological Institute.

%prep
rm -rf $RPM_BUILD_ROOT

%setup -q -n %{NAME}

%build
make %{_smp_mflags}

%install
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig


%files
%defattr(0775,root,root,0775)
%{_libdir}/lib%{NAME}.so

%package -n %{NAME}-devel
Summary: Road surface weather model development files
Provides: %{NAME}-devel
Requires:%{NAME}

%description -n %{NAME}-devel
Road surface weather model development files

%files -n %{NAME}-devel
%defattr(0664,root,root,0775)
%{_includedir}/roadsurf/roadsurf.mod
%{_includedir}/roadsurf/roadsurfvariables.mod
%{_includedir}/roadsurf/Constants.h

%changelog
* Thu Jul  3 2025 Virve Karsisto <virve.karisto@fmi.fi> 25.7.3-1.fmi
- -h

* Tue Feb 11 2025 Mika Heiskanen <mika.heiskanen@fmi.fi> - 1.4.1-1.fmi
- bug fix to prevent water storage increasing by double
- Bug fix to prevent melting huge amounts at once

* Mon Aug 19 2024 Andris Pavenis <andris.pavenis@fmi.fi> 1.3.1-1.fmi
- Changed output array initialization from -99 to -9999 (update from Virve Karsisto <virve.karsisto@fmi.fi>)

* Fri Jul 12 2024 Andris Pavēnis <andris.pavenis@fmi.fi> 1.3-1.fmi
- Replace many boost library types with C++ standard library ones

* Wed Nov 15 2023 Virve Karsisto <tierut@dev.road.fmi.fi>
- Bug fix to secondary ice storage wear

* Wed Sep 20 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 1.1-1.fmi
- Number of ground layers can now be changed from configuration file
- Fixed error in getTempAtDepth arguments
- Fixed error in melting

* Mon Jun 12 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 1.0-1.fmi
- Release version

* Wed May 10 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.3-1.fmi
- Improved class and method names

* Tue May  2 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.2-1.fmi
- Removed obsolete iterative heat balance code

* Tue May  2 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.1-1.fmi
- Initial open source version
