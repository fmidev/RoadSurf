%define BINNAME roadrunner
%define RPMNAME smartmet-%{BINNAME}
Summary: Road weather model
Name: %{RPMNAME}
Version: 23.5.2
Release: 1%{?dist}.fmi
License: FMI
Group: Development/Tools
URL: https://github.com/fmidev/smartmet-roadrunner
Source0: %{name}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot-%(%{__id_u} -roadn)

%if 0%{?rhel} && 0%{rhel} < 9
%define smartmet_boost boost169
%else
%define smartmet_boost boost
%endif

%define smartmet_fmt_min 8.1.1
%define smartmet_fmt_max 8.2.0

BuildRequires: roadsurf-devel
BuildRequires: rpm-build
BuildRequires: gcc-c++
BuildRequires: make
BuildRequires: %{smartmet_boost}-devel
BuildRequires: smartmet-library-newbase-devel >= 22.11.14
BuildRequires: smartmet-library-macgyver-devel >= 22.12.16
BuildRequires: smartmet-library-locus-devel >= 22.12.16
BuildRequires: libstx-exparser-devel >= 22.2.14
BuildRequires: gdal34-devel
BuildRequires: cpr
BuildRequires: jsoncpp-devel
BuildRequires: libcurl-devel
BuildRequires: fmt-devel >= %{smartmet_fmt_min}, fmt-devel < %{smartmet_fmt_max}
Requires: smartmet-library-newbase >= 22.11.14
Requires: smartmet-library-macgyver >= 22.12.16
Requires: smartmet-library-locus >= 22.12.16
Requires: libstx-exparser >= 22.2.14
Requires: roadsurf
Requires: %{smartmet_boost}-regex
Requires: %{smartmet_boost}-iostreams
Requires: %{smartmet_boost}-filesystem
Requires: %{smartmet_boost}-program-options
Requires: %{smartmet_boost}-system
Requires: jsoncpp
Requires: libcurl
Requires: gdal34-libs
Requires: fmt >= %{smartmet_fmt_min}, fmt < %{smartmet_fmt_max}
Provides: roadrunner


%description
Road weather modelBuild
Requires: libstx-exparser-devel >= 22.2.14

%prep
rm -rf $RPM_BUILD_ROOT

%setup -q -n %{RPMNAME}

%build
make %{_smp_mflags}

%install
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(0775,root,root,-)
%{_bindir}/roadrunner

%changelog
* Tue May  2 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 23.5.2-1.fmi
- Use open source roadsurf library
* Thu Jan 12 2023 Mika Heiskanen <mika.heiskanen@fmi.fi> - 23.1.12-1.fmi
- Optimized snowfall calculations for speed
* Wed Jan 11 2023 Mika Heiskanen <mika.heiskanen@fmi.fi> - 23.1.11-1.fmi
- Fixed array overflow issue
- Thread safety fix
* Tue Jan 10 2023 Mika Heiskanen <mika.heiskanen@fmi.fi> - 23.1.10-1.fmi
- Fixes to ploughing algorithms
* Mon Jan  9 2023 Mika Heiskanen <mika.heiskanen@fmi.fi> - 23.1.9-1.fmi
- Improvements to sky view runs
* Wed Dec 21 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.12.21-1.fmi
- Thread safety fixes
- Improved surface layer temperature initialization from road weather observations
- Fixed issues reported by CodeChecker
* Fri Dec 16 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.12.16-1.fmi
- Repackaged since PostgreSQLConnection ABI changed
* Tue Oct 25 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.10.25-1.fmi
- Changed frostsum to use observed tsurf when available
* Mon Oct 17 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.10.17-1.fmi
- Fixes to FrostSum calculation
- Speed optimizations
* Wed Oct  5 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.10.5-1.fmi
- Fixed FrostSum calculation
* Mon Sep 26 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.9.26-2.fmi
- Less verbose output
* Mon Sep 26 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.9.26-1.fmi
- Fixed missing value handling of frost sum calculations
- Less verbose output
* Mon Jun 20 2022 Andris Pavēnis <andris.pavenis@fmi.fi> 22.6.20-1.fmi
- Add support for RHEL9, upgrade libpqxx to 7.7.0 (rhel8+) and fmt to 8.1.1
* Wed Jun  8 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.6.8-1.fmi
- Added frost sum calculation
- Added optional sky view input parameters
* Tue May 24 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.5.24-1.fmi
- Repackaged due to NFmiArea ABI changes
* Fri May 20 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.5.20-1.fmi
- Repackaged due to ABI changes to newbase LatLon methods
* Thu Feb 17 2022 Andris Pavēnis <andris.pavenis@fmi.fi> 22.2.17-1.fmi
- Support both RHEL7 and RHEL8
* Mon Feb 14 2022 Andris Pavēnis <andris.pavenis@fmi.fi> 22.2.14-1.fmi
- Use new build of stx-exparser and update dependencies
* Mon Sep 20 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.9.20-1.fmi
- Added RoadNotification
* Fri Sep 17 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.9.17-1.fmi
- Added traffic index reasoning
- Added traffic reasoning logic to the model and replaced RoadReasoning output with it.
- Added output variable RoadConditionSeverity that is now Traffic value output instead
* Tue Jun 29 2021 Pertti Kinnia <pertti.kinnia@fmi.fi> - 21.6.29-1.fmi
- Maintenance mode updates
- Bug fixes and update to salting conditions
- Fixed bug with freesurf storages and forced everything to melt when salt on road
* Fri Jun 18 2021 Pertti Kinnia <pertti.kinnia@fmi.fi> - 21.6.18-1.fmi
- Added ability to set the depth from which the output surface temperature is interpolated in config file
* Thu May  6 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.5.6-1.fmi
- Repackaged due to ABI changes in NFmiAzimuthalAre
* Tue Mar  2 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.3.2-1.fmi
- Repackaged with latest base libraries
* Fri Feb 19 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.2.19-1.fmi
- Improved radiation handling
* Mon Jan 18 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.1.18-1.fmi
- Fixed himan precipitation interpretations
* Tue Dec 15 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.12.15-1.fmi
- Upgrade to pgdg12
* Wed Oct 28 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.10.28-1.fmi
- Upgrade to fmt 7.1
* Fri Aug 21 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.8.21-1.fmi
- Upgrade to fmt 6.2
* Mon Jun  8 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.6.8-1.fmi
- Bug fixes
- Upgraded libpqxx
* Mon Apr 20 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.20-1.fmi
- Read the coupling period length from the configuration file
* Sat Apr 18 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.18-1.fmi
- Upgraded to Boost 1.69
* Thu Apr 16 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.16-1.fmi
- Added querydata masking to limit the number of simulation points for example at water points
* Thu Dec  5 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.12.5-1.fmi
- Optimized Fortran pow calls for a 50% speed increase
* Wed Dec  4 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.12.4-2.fmi
- Reduced NFmiMetTime constructions for a 25% speed increase
* Wed Dec  4 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.12.4-1.fmi
- Added checks for the case when querydata times do not cover the simulation interval at any point
* Wed Nov 20 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.11.20-1.fmi
- Repackaged due to newbase API changes
* Thu Nov  7 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.11.7-1.fmi
- Fixed precipitation rate parameter name typo
* Thu Oct 31 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.31-1.fmi
- Rebuilt due to newbase API/ABI changes
* Fri Oct 25 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.25-1.fmi
- Removed debugging prints
* Tue Oct 22 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.22-1.fmi
- Fixed querydata time indexing
* Tue Oct 15 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.15-1.fmi
- Fixed nearest station distance conditions
* Fri Oct 11 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.11-1.fmi
- Added possibility to list wanted fmisid:s for SmartMet Server requests
* Fri Sep 27 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.9.27-1.fmi
- Repackaged due to ABI changes in SmartMet libraries
* Mon Sep 23 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.9.23-1.fmi
- Bug fix to relaxation
- Calculate humidity from tdew
* Wed Aug 28 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.8.28-1.fmi
- Repackaged since Locus Location ABI changed
* Tue Jul 16 2019 Pertti Kinnia <pertti.kinnia@fmi.fi> - 19.7.16-1.fmi
- Bug fix, new release version
* Tue Jul  2 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.7.2-1.fmi
- Fixes to coupling
* Tue May 21 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.5.21-2.fmi
- Fixed PrecipitationRate to work
* Tue May 21 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.5.21-1.fmi
- Fixed thread_local storage for querydata input
* Thu Mar 21 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.3.21-1.fmi
- Added option --infile (-i) to define an ASCII data source from the command line
- Option -o now works for single point research runs too, default is still stdout
- JSON config input.points can now list one or more ASCII data source filenames
* Wed Jan 30 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.1.30-1.fmi
- Fixed option -o to work in parallel mode too
* Fri Jan 25 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.1.25-1.fmi
- Added reading of maintenance operations
- Added option --outfile (-o) for setting the output filename
* Thu Dec 13 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.12.13-1.fmi
- First release of the joint Fortran/C++ version
