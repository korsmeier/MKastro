
__author__='M. Korsmeier - (korsmeier@physik.rwth-aachen.de)'
__version__='1.0'


  ###############
 #  Constants  #
###############

# speed of light
c = {'m/s':299792458., 'kpc/My':306.3915}

# nuclei masses
m_proton  = {'GeV':0.93827, 'MeV':938.27}
m_neutron = {'GeV':0.93827, 'MeV':938.27}

Z_dictionary = {1:'proton',
                2:'helium',
                3:'lithium',
                4:'beryllium',
                5:'boron',
                6:'carbon',
                7:'nitrogen',
                8:'oxigen',
                9:'flour',
                10:'neon',
                11:'natrium',
                12:'magnesium',
                13:'aluminium',
                14:'silicon',
                -1:'antiproton',
                -2:'antihelium'
            }

  #####################
 #  Transformations  #
#####################

length = {
    'cm->kpc': 3.24078e-22,
     'm->kpc': 3.24078e-20,
    'km->kpc': 3.24078e-17,

    'cm->pc' : 3.24078e-19,
     'm->pc' : 3.24078e-17,
    'km->pc' : 3.24078e-14,

    'kpc->cm': 3.086e+21,
    'kpc->m' : 3.086e+19,
    'kpc->km': 3.086e+16,

    'pc->cm' : 3.086e+18,
    'pc->m'  : 3.086e+16,
    'pc->km' : 3.086e+13,

    'cm->m'  : 1e-2,
    'cm->km' : 1e-5,
        
    'km->cm' : 1e5,
    'km->m'  : 1e3,

    'm->cm'  : 1e2,
    'm->km'  : 1e-3
}

time  = {
    's->min': 60.,
    's->h'  : 3600.,
    's->y'  : 3.17098e-8,
    's->My' : 3.17098e-14,

    'y->s'  : 3.154e+7,
    'My->s' : 3.154e+13
}

energy  = {
    'GeV->TeV': 1e-3,
    'TeV->GeV': 1e3,
    'GeV->MeV': 1e3,
    'MeV->GeV': 1e-3

}

diffusion_constant = {
    'cm^2/s->kpc^2/My': length['cm->kpc']**2/time['s->My'],
    'kpc^2/My->cm^2/s': length['kpc->cm']**2/time['My->s']
}

flux = {
    'MeV^-1cm^-2->GeV^-1m^-2': 1./length['cm->m']**2/energy['MeV->GeV'],
    'GeV^-1m^-2->MeV^-1cm^-2': 1./length['m->cm']**2/energy['GeV->MeV'],
    'kpc^-2My^-1->m^-2s^-1': 1./length['kpc->m']**2/time['My->s'],
    'm^-2s^-1->kpc^-2My^-1': 1./length['m->kpc']**2/time['s->My']
}

XS ={
    'b->m^2'   : 1e-28,
    'mb->m^2'  : 1e-31,
    'b->cm^2'  : 1e-24,
    'mb->cm^2' : 1e-27,

    'm^2->b'   : 1e28,
    'm^2->mb'  : 1e31,
    'cm^2->b'  : 1e24,
    'cm^2->mb' : 1e27,
}
