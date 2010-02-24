//  Thomas-Fermi kinetic energy functional

template<class num>
static num thomasfermi_kinetic(const densvars<num> &d)
{
    using xc_constants::CF;

    return CF*pow(d.n, 5.0/3.0);
}

