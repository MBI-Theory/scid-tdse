/^#/{next}
{ t  = $(ct) ; f = $(cf) ; g = $(cg) ;
  if (t>ot) {
    gn = (f-of)/(t-ot) ;
    printf "%20.12g %20.12g %20.12g\n", t, gn, 0.5*(g+og) ;
    }
  of = f ; ot = t ; og = g ;
  }
