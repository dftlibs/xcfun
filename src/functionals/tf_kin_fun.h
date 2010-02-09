//  Thomas-Fermi kinetic energy functional

template<class num>
static num thomasfermi_kinetic(const num &rhoa,
			       const num &rhob)
{
    using xc_constants::CF;

    return CF*pow(rhoa+rhob, 5.0/3.0);
}

