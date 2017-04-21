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

  //Analysis::CombinedAnalysis analysis(mcp, pot, spc);
  Histogram<double> cmdist(1);
  spc.load("state");                                  // Load start configuration, if any

  Group agg = *spc.findMolecules("fibril")[0];

  agg.translate( spc, -Geometry::massCenter( spc.geo, spc.p, agg) );
  agg.accept( spc );

  for (auto i : agg )
    if (atom[spc.p[i].id].name == "CTR")
      cout << atom[spc.p[i].id].name << endl;
  
  FormatPQR::save("init.pqr", spc.p, spc.geo.len);

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      mv.move();
       // double dist = gouychapman->expot.surfDist( Geometry::massCenter( spc.geo, spc.p, agg) );
       //  cmdist( dist )++;
    } // end of micro loop

    cout << loop.timing();  // print timing and ETA information

    spc.save("state");               
    //analysis.sample();
    //cmdist.save("cmdist.dat");

  } // end of macro loop
double min = mcp["energy"]["penalty"]["xyz"]["lo1"];
penalty->save("pf_", 1, {min+5, min+15} );
cout << loop.info() + mv.info() + spc.info() + penalty->info();
}
