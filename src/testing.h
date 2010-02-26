
/* Run a second order calculation of functional fun_id,
   compare to the reference values. Return the number
   of values in error. Note that the density
   values should be in order

   density: na nb saa sbb sab
   reference: e de/dna de/dnb de/dsaa de/dsbb de/dsab
              na.na na.nb na.saa na.sbb na.sab
	      nb.nb nb.saa nb.sbb nb.sab
	      saa.saa saa.sbb saa.sab 
	      sbb.sbb sbb.sab sab.sab


   calc - ref + x*eref*ref = 0 => (ref - calc) = x*ref*eref
   x is in [-1:1] if fabs(ref - calc) < fabs(ref*eref)
   
*/
template<int fun_id>
int standard_abgga_test(const double density[5], 
			const double reference[21],
			double ref_rel_err)
{
  densvars<taylor<double,5,2> > d;
  setdens<double,XC_ABGGA,2>(d,density);
  taylor<double,5,2> res = functional<fun_id>::energy(d);
  res.deriv_facs();
  int nerr = 0;

  for (int i=0;i<res.size;i++)
    if (fabs(res[i] - reference[i]) > fabs(reference[i]*ref_rel_err))
      nerr++;

#ifndef NO_STDCXX
  cout << scientific;
  cout.precision(16);
  if (nerr > 0)
    {
      cout << "Error detected in functional " <<
	functional<fun_id>::get_name() << endl;
      cout << "Abs.Error \tComputed              Reference" << endl;
      for (int i=0;i<res.size;i++)
	{
	  cout.precision(3);
	  cout << fabs(res[i]-reference[i]);
	  cout.precision(16);
	  cout << "\t " << res[i] << " " << reference[i];
	  if (fabs(res[i] - reference[i]) > fabs(reference[i]*ref_rel_err))
	    cout << " *" << endl;
	  else
	    cout << endl;
	}
    }
#endif
  return nerr;
}
