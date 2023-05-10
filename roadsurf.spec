%define NAME roadsurf
Summary: Road surface weather model
Name: %{NAME}
Version: 0.1.3
Release: 1%{?dist}.fmi
License: GPL
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
* Wed May 10 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.3-1.fmi
- Improved class and method names

* Tue May  2 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.2-1.fmi
- Removed obsolete iterative heat balance code

* Tue May  2 2023 Mika Heiskanen <mheiskan@rhel8.dev.fmi.fi> - 0.1.1-1.fmi
- Initial open source version
