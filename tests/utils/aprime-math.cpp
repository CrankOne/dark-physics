# include <iostream>
# include <fstream>
# include <cmath>
# include <set>
# include <cstring>
# include <unistd.h>
# include <cassert>
# include <tuple>

# include "dphmc-test-aprime-ww.hh"

namespace dphmc {
namespace test {

/// Generate reference curves for \chi/Z^2 on Fig. 10 of \cite{Bjorken}.
void
bjorken_chi_ref_curves( const std::string & outputDir ) {
    struct {
        std::string label;
        double values[3];
    } ps[] = {
        { "W-74-183-.2",    {  0.2, 183.84, 74 } },
        { "W-74-183-1",     {    1, 183.84, 74 } },
        { "W-74-183-6",     {    6, 183.84, 74 } },
        { "Al-13-26-6",     {    6,  26.98, 13 } },
        { "Pb-82-207-80",   {   80, 207.2,  82 } },
        { "", { -1 } }
    };
    for( auto p = ps; !p->label.empty(); ++p ) {
        std::ofstream os = std::ofstream( outputDir + p->label + ".dat" );
        std::cout << "xxx" << std::endl;  // XXX
        os << "# " << p->label << std::endl;
        double tstVal;
        for( double am = 1e-3; am <= 6. && am <= p->values[0] ; am *= 1.5 ) {
            try {
                tstVal = chi_bjorken( &am, p->values );
            } catch( std::runtime_error &e ) {
                tstVal = std::nan("");
            }
            os << am << " " << tstVal << std::endl;
            if( std::isnan(tstVal) ) {
                break;
            }
        }
        os.close();
    }
}

void
aprime_dcs_surf( const struct dphmc_APrimePhysParameters & ps
               , struct dphmc_IterativeQAGSParameters & it
               , FILE * os
               , size_t nPointsX, size_t nPointsTheta
               , FILE * mxOs ) {
    fprintf( os, "# A' WW-approach configuration: \n# " );
    dphmc_print_aprime_ws_parameters( os, &ps );
    fprintf( os, "\n# number of points x, theta: %zu, %zu.\n"
           , nPointsX, nPointsTheta );
    double xRange[2], thetaRange[2];
    auto points = aprime_dcs_ww( ps, it, nPointsX, nPointsTheta
                                            , xRange, thetaRange );
    fprintf( os, "# x range: %e, %e\n", xRange[0], xRange[1] );
    fprintf( os, "# theta range: %e, %e\n", thetaRange[0], thetaRange[1] );
    fprintf( os, "# columns: x, theta, d \\sigma, d^2 \\sigma/d \\theta: s/a, -\"-: num, -\"-: num.err\n" );
    size_t xCounter = 0;
    for( auto pt : points ) {
        // Write x, theta and diff sigma values
        fprintf( os, "%e\t%e\t%e\t", pt.x, pt.theta, pt.dCS, pt.saDDTheta );
        // Write semi-analytic and numerical estimation of diff sigma by theta
        fprintf( os, "%e\t%e\t%e\n", pt.saDDTheta, pt.numDDTheta, pt.numDDThetaErr );
        if( ! ((++xCounter) % nPointsTheta) )
            fprintf( os, "\n" );
    }
    if( mxOs ) {
        auto pts = aprime_dcs_max( ps, it, nPointsX );
        for( auto pt : pts ) {
            fprintf( mxOs, "%e\t%e\t%e\n", std::get<0>(pt), std::get<1>(pt), std::get<2>(pt) );
        }
    }
}

int
aprime_full_cs( const struct dphmc_APrimePhysParameters & ps
              , struct dphmc_IterativeQAGSParameters & it
              , size_t nCalls
              , double (*  elastic_f)( double, struct dphmc_APrimePhysParameters * )
              , double (*inelastic_f)( double, struct dphmc_APrimePhysParameters * )
              , int flags
              ) {
    fputs( "# ", stdout );
    dphmc_print_aprime_ws_parameters( stdout, &ps );
    printf("\n# VEGAS integration; N=%zu\n", nCalls );
    # if 0
    void * ws = NULL;
    int rc = dphmc_init_aprime_cs_workspace( &ps, &ws );
    if( rc ) {
        fprintf( stderr
               , "Failed to initialize A' cross-section workspace: %d.\n"
               , rc );
        dphmc_free_aprime_cs_workspace( ws );
        return -1;
    }
    #endif
    aprime::WWCaches * wwCaches = dphmc_aprime_new( &ps
                             , &it
                             , elastic_f
                             , inelastic_f
                             , flags
                             );

    assert(wwCaches);
    double relErrPtr, absErr
         , chi2
         , result1 = dphmc_aprime_ww_full_numeric_1( &it, wwCaches, &relErrPtr, &absErr)
         , result2 = dphmc_aprime_ww_full_numeric_2( wwCaches
                                                   , nCalls
                                                   , &absErr
                                                   , &chi2 );
    std::cout << "estimation = "
              << dphmc_aprime_ww_fast_integral_estimation(wwCaches)
              << std::endl;
    //dphmc_free_aprime_cs_workspace( ws );
    dphmc_aprime_delete( wwCaches );
    std::cout << "full CS1 = " << result1 << std::endl;
    std::cout << "CS1 rel err. = " << relErrPtr << std::endl;
    std::cout << "full CS2 = " << result2 << std::endl;
    std::cout << "abs. error = " << absErr << std::endl;
    std::cout << "chi2 = " << chi2 << std::endl;
    return 0;
}

}  // namespace ::dphmc::test
}  // namespace ddphmc

static void
_print_usage( const char * appName, std::ostream & os ) {
    os << "This application assists for the numerical tests of the \"dark"
        " physics\" Monte-Carlo package on A' production. Usage:" << std::endl
       << "    $ " << appName << " [-h|--help] -- prints this message to stdout."
       << std::endl;
    os << "    $ " << appName << " chi <output-dir> -- writes reference points"
        " on chi curves (demon purposes) at given dir"
       << std::endl;
    os << "    $ " << appName << " acs <output-file> <n-points-X> <n-points-theta> "
        " [-Z <materialZ=74>]"
        " [-A <materialA=183.84>]"
        " [-m <A'-mass-GeV=0.05>] "
        " [-b <beam-energy-GeV=6>]"
        " [-e <epsilon=1e-4>]"
        " [-a <GK-absolute-error=1e-12>]"
        " [-r <GK-relative-error=1e-12>]"
        " [-i <GK-increment=1.1>]"
        " [-l <GK-limit=1000>]"
        " [-n <GK-nnodes=1000>]"
        " [-X <max-output-file.dat>]"
        " [-f <flag-name>]" << std::endl <<
        " Application writes reference points"
        " on differential cross-section surface for given parameters of"
        " calculation to the output dir. The -A and -Z set target atomic"
        " numbers, the -b sets beam energy in GeV, -e sets the mixing constant,"
        " and -m sets the A' hypothetical mass in GeV. Note that procedure uses"
        " GSL's Gauss-Kronrod integration method, so GK-prefixed parameters are"
        " related to it. See \"Numerical Integration\" section of `nfo gsl`"
        " for detailed reference. If -X option provided, its argument is"
        " interpreted as output filename for theta_max value (at some x)."
        " A' WW approximation cross section procedure may be additionally"
        " steered by flags: \"enable-chi-theta-dep\","
        " \"disable-theta-cut-off\", \"disable-phys-parameter-checks\"."
        " Corresponding documentation for flags may be found at the online doc"
        " for function dphmc_aprime_new()."
       << std::endl;
    os << "    $ " << appName << " int ..."
       << std::endl;
}

static int
_configure_aprime_pars_from_command_line( int argc, char * const argv[]
                                        , struct dphmc_APrimePhysParameters & ps
                                        , struct dphmc_IterativeQAGSParameters & it
                                        , size_t * nCalls
                                        , char ** maxOutFileNamePtr
                                        , int & flags ) {
    int c;
    opterr = 0;
    flags = 0x0;
    while ((c = getopt(argc, argv, "Z:A:b:m:e:a:r:i:l:n:N:X:f:")) != -1)
        switch (c) {
            case 'Z': ps.Z = (uint16_t) atoi( optarg ); break;
            case 'A': ps.A = atof( optarg );            break;
            case 'm': ps.massA_GeV = atof( optarg );    break;
            case 'b': ps.EBeam_GeV = atof( optarg );    break;
            case 'e': ps.epsilon = atof( optarg );      break;
            case 'a': it.epsabs = atof( optarg );       break;
            case 'r': it.epsrel = atof( optarg );       break;
            case 'i': it.epsrelIncFt = atof( optarg );  break;
            case 'l': it.limit = atoi( optarg );        break;
            case 'n': it.nnodes = atoi( optarg );       break;
            // A' WW calc flags -- see doc on dphmc_aprime_new() flag's
            // argument
            case 'f': {
                if( !strcmp( optarg, "enable-chi-theta-dep" ) ) {
                    flags |= 0x1;
                } else if( !strcmp( optarg, "disable-theta-cut-off" ) ) {
                    flags |= 0x2;
                } else if( !strcmp( optarg, "disable-phys-parameter-checks" ) ) {
                    flags |= 0x4;
                } else {
                    std::cerr << "Unrecognized A' WW calculation flag \""
                              << optarg
                              << "\"." << std::endl;
                    return 1;
                }
            } break;
            case 'N': {
                if( !nCalls ) {
                    std::cerr << "-N is valid for `int' procedure only."
                              << std::endl;
                    return 1;
                }
                *nCalls = atoi( optarg );
            } break;
            case 'X' : {
                if( !maxOutFileNamePtr ) {
                    std::cerr << "-X is valid for `acs' procedure only."
                              << std::endl;
                    return 1;
                }
                *maxOutFileNamePtr = optarg;
            } break;
            case '?':
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                return 1;
    }
    return 0;
}

int
main(int argc, char * const argv[] ) {
    std::set<std::string> args( argv + 1, argv + argc );
    if( 0 == args.size()
     || args.end() != args.find("-h")
     || args.end() != args.find("--help") ) {
        _print_usage( argv[0], std::cout );
        return 0;
    }
    if( argc > 1 ) {
        if( !strcmp( argv[1], "chi" ) ) {
            if( argc != 3 ) {
                std::cerr << "Error: argument <output-dir> is missing."
                          << std::endl;
                _print_usage( argv[0], std::cerr );
                return -1;
            }
            dphmc::test::bjorken_chi_ref_curves( argv[2] );
            return EXIT_SUCCESS;
        }
        bool acs;
        if( (acs = !strcmp( argv[1], "acs" ))
         || !strcmp( argv[1], "int" ) ) {
            struct dphmc_APrimePhysParameters ps {
                /* Z ............... */ (uint16_t) 74,
                /* A ............... */ 183.84,
                /* mass A', GeV .... */ 0.05,
                /* E beam, GeV ..... */ 6,
                /* mixing factor ... */ 1e-4,
                /* re-norm factor .. */ 1,
                /* theta cut-off c . */ 100  // TODO: check this
            };
            struct dphmc_IterativeQAGSParameters it {
                /* epsabs .......... */ 1e-12,
                /* epsrel .......... */ 1e-12,
                /* epsrelIncFt ..... */ 1.1,
                /* limit ........... */ (size_t) 1e3,
                /* nnodes .......... */ (size_t) 1e3
            };
            int flags;
            if( acs && argc < 5 ) {
                std::cerr << "Error: output file path and # of &theta points"
                    " must be specified as first positional arguments after"
                    " procedure name \"acs\"."
                    << std::endl;
                _print_usage( argv[0], std::cerr );
                return -1;
            }
            size_t nPointsX, nPointsTheta;
            size_t nCalls = 1e5;
            int cc;
            char * maxOutFileName = NULL;
            if(acs) {
                nPointsX = atoi( argv[3] );
                nPointsTheta = atoi( argv[4] );
                cc = _configure_aprime_pars_from_command_line(
                          argc - 4, argv + 4
                        , ps
                        , it
                        , &nCalls
                        , &maxOutFileName
                        , flags );
            } else {
                cc = _configure_aprime_pars_from_command_line(
                          argc - 1, argv + 1
                        , ps
                        , it
                        , &nCalls
                        , NULL
                        , flags );
            }
            if( cc ) {
                std::cerr << "Error: exit due to previous configuration error."
                          << std::endl;
                return EXIT_FAILURE;
            }
            if( !acs ) {
                #if 0
                throw std::runtime_error( "TODO: full CS with (in)elastic" );
                #else
                int rc = ::dphmc::test::aprime_full_cs( ps, it, nCalls
                                                      , dphmc_aprime_form_factor_elastic
                                                      , dphmc_aprime_form_factor_inelastic
                                                      , flags
                                                      );
                if( rc ) return EXIT_FAILURE;
                #endif
                return EXIT_SUCCESS;
            } else {
                FILE * osF = fopen( argv[2], "w" )
                   , * mxOsF = NULL;
                if( maxOutFileName ) {
                    mxOsF = fopen( maxOutFileName, "w" );
                }
                ::dphmc::test::aprime_dcs_surf( ps, it
                                              , osF, nPointsX
                                              , nPointsTheta, mxOsF );
                fclose( osF );
                if( mxOsF ) {
                    fclose( mxOsF );
                }
            }
            return EXIT_SUCCESS;
        }
        std::cerr << "Argument #1 \"" << argv[1] << "\" is not recognized."
                  << std::endl;
        _print_usage( argv[0], std::cerr );
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
