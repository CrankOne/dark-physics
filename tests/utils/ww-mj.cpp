

void
usage_info( std::ostream & os, const std::string & appName ) {
    os << "Usage:" << std::endl
       << "    $ " << appName << " [1|2] <outFile>" << std::endl
       << "Application samples the WW cross section values for certain physics"
       << " and writes the resulting 100-bin histogram into a file."
       << std::endl;
}

int
main(int argc, char * argv[]) {
    return EXIT_SUCCESS;
}
