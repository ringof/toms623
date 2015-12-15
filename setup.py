from numpy.distutils.core  import setup, Extension

# interface for Renka's algorithm 623 fortran code
ext = Extension(name  = '_toms623',
                sources       = ['_toms623.pyf','_toms623.f90'])

if __name__ == "__main__":
    setup(name = 'toms623',
          version           = "0.0.1",
          description       = "Python interface to TOMS 623 fortran code",
          author            = "Jeff Whitaker",
          author_email      = "jeffrey.s.whitaker@noaa.gov",
          ext_modules       = [ext],
          packages          = ['toms623'],
          )
