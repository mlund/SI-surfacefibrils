#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboidslit> Tspace;

int main() {
  InputMap mcp("input.json");         // Open input parameter file
  MCLoop loop(mcp);                  // handle mc loops
  Tspace spc(mcp); 

  auto pot = Energy::ExternalPotential<Tspace,Potential::GouyChapman<double,true> >(mcp)
    + Energy::PenaltyEnergy<Tspace>(mcp, spc);
  auto gouychapman = std::get<0>(pot.tuple());
  auto penalty = std::get<1>(pot.tuple());

  gouychapman->expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

  // MC moves
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  penalty->load("pf_");

  double lo1 = mcp["energy"]["penalty"]["xyz"]["lo1"];
  double hi1 = mcp["energy"]["penalty"]["xyz"]["hi1"];
  std::vector<double> bw = { 1, 1 };
  std::vector<double> lo = { lo1, 0 };
  std::vector<double> hi = { hi1, 90 };
  Table<double> hist2d( bw, lo, hi );
  cout << lo[0] << " " << hi[0] << endl;
  spc.load("state");                                  // Load start configuration, if any

  Group agg = *spc.findMolecules("fibril")[0];

  agg.translate( spc, -Geometry::massCenter( spc.geo, spc.p, agg) );
  agg.accept( spc );

  for (auto i : agg )
    if (atom[spc.p[i].id].name == "CTR")
      cout << atom[spc.p[i].id].name << endl;
  
  FormatPQR::save("init.pqr", spc.p, spc.geo.len);

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
  Point dir;
  dir << (mcp["energy"]["penalty"]["xyz"]["dir"] | std::string("0 0 1"));
  while ( loop[0] ) {  // Markov chain
    while ( loop[1] ) {
      mv.move();
      if ( slump() > .5) {
        auto cm = Geometry::massCenter(spc.geo, spc.p, agg);
        double dist = cm.z(); 
        auto S = Geometry::gyration(spc.geo, spc.p, agg, cm);
        Eigen::EigenSolver<Eigen::Matrix3d> es(S);
        double max = std::numeric_limits<double>::min();
        int j = 0;
        for ( int i = 0; i < es.eigenvalues().size(); ++i )
        {
            double value = es.eigenvalues()[i].real();
            if ( value > max )
            {
                max = value;
                j = i;
            }
        }
        Point eig = es.eigenvectors().col(j).real();
        double angle = acos( eig.dot(dir) ) * 180. / pc::pi;
        if ( angle > 90 )
            angle = 180 - angle;
        std::vector<double> v = { dist, angle };
        hist2d.to_index(v);
        hist2d[ v ] += 1;
      }
    } // end of micro loop

    cout << loop.timing();  // print timing and ETA information

    spc.save("state");               
    hist2d.save("hist2d.dat");

  } // end of macro loop
cout << loop.info() + mv.info() + spc.info() + penalty->info();
}
