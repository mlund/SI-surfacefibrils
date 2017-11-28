#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cylinder> Tspace;

int main(int argc, char** argv) {

  InputMap mcp("input.json");

  Tspace spc(mcp);
  
  Group agg = *spc.findMolecules("fibril")[0];
  Group sod = *spc.findMolecules("sodium")[0];
  
  // mass center at the origin
  agg.translate( spc, -Geometry::massCenter( spc.geo, spc.p, agg) );
  agg.accept( spc );
  
  // align aggregate principal axis to cylinder axis
  auto cm = Geometry::massCenter(spc.geo, spc.p, agg);
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
  Point dir = {0, 0, 1};
  double angle = acos( eig.dot( dir ) );
  cout << "angle between aggregate and z-axis: " << angle * 180 / pc::pi << endl;

  auto endpoint = eig.cross( dir );
  agg.rotate( spc, endpoint, angle );
  agg.accept( spc );

  // calculate the angle between aggregate and z-axis after rotation
  cm = Geometry::massCenter(spc.geo, spc.p, agg);
  S = Geometry::gyration(spc.geo, spc.p, agg, cm);
  Eigen::EigenSolver<Eigen::Matrix3d> es2(S);
  max = std::numeric_limits<double>::min();
  j = 0;
  for ( int i = 0; i < es2.eigenvalues().size(); ++i )
  {
      double value = es2.eigenvalues()[i].real();
      if ( value > max )
      {
          max = value;
          j = i;
      }
  }
  eig = es2.eigenvectors().col(j).real();
  dir = {0, 0, 1};
  angle = acos( eig.dot( dir ) );
  cout << "angle between aggregate and z-axis: " << angle * 180 / pc::pi << endl;

  // displace sodium ions from the fibril region
  for ( auto i : sod ) {
    if ( sqrt(spc.p[i].x()*spc.p[i].x()+spc.p[i].y()*spc.p[i].y()) < 55 ) {
      double x = spc.p[i].x();
      double diffx = (x > 0) ? 55-x : -55-x;
      double y = spc.p[i].y();
      double diffy = (y > 0) ? 55-y : -55-y;
      Point t = {diffx, diffy, 0};
      spc.trial[i].translate(spc.geo, t);
      spc.p[i] = spc.trial[i];
      sod.accept( spc );
    }
  }

  // save initial configuration
  FormatPQR::save("init.pqr", spc.p, spc.geo.len);

  double len = mcp["system"]["geometry"]["length"];
  double rad = mcp["system"]["geometry"]["radius"];
  int macro = mcp["system"]["mcloop"]["macro"];
  int micro = mcp["system"]["mcloop"]["micro"];
  int tot = macro * micro;

  cout << "cylinder length: " << len << ", radius: " << rad << endl;

  //std::vector<double> bw = { 1, 1 };
  //std::vector<double> lo = { -len/2., 0 };
  //std::vector<double> hi = { len/2., rad };
  //Table<double> histFib( bw, lo, hi );
  std::vector<double> bw = { .1, 1 };
  std::vector<double> lo = { 0, 0 };
  std::vector<double> hi = { rad, 0 };
  Table<double> histChar( bw, lo, hi );
  Table<double> histCat( bw, lo, hi );
  Table<double> histCharAgg( bw, lo, hi );
  for ( auto i : agg ) {
    double radial = sqrt( spc.p[i].x()*spc.p[i].x() + spc.p[i].y()*spc.p[i].y() ); 
    std::vector<double> v = { radial };
    histCharAgg.to_index(v);
    histCharAgg[ v ] += atom[spc.p[i].id].charge;
  }
  histCharAgg.save("histCharAgg.dat");

  double z_min = len/2.;
  double z_max = -len/2.;
  for ( auto i : agg ) 
  {
    //double radial = sqrt( spc.p[i].x()*spc.p[i].x() + spc.p[i].y()*spc.p[i].y() ); 
    double axial = spc.p[i].z();
    if (axial > z_max)
      z_max = axial;
    if (axial < z_min)
      z_min = axial;
    //std::vector<double> v = { axial, radial };
    //histFib.to_index(v);
    //histFib[ v ] += 1;
  }
  //histFib.save("histFib.dat");
  cout << "z_min: " << z_min << ", z_max: " << z_max << endl;

  auto pot = Energy::Nonbonded<Tspace,CoulombWCA>(mcp);
  Analysis::CombinedAnalysis analysis(mcp,pot,spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  int sweeps = 0;
  int sample = 0;
  while ( loop[0] ) {  // Markov chain
    while ( loop[1] ) {
      sweeps++;
      mv.move();
      analysis.sample();
      if ( sweeps > 5e4 && slump() > .5) {
        sample++;
        for ( auto i : sod ) {
          if ( spc.p[i].z() > z_min && spc.p[i].z() < z_max ) {
            double radial = sqrt( spc.p[i].x()*spc.p[i].x() + spc.p[i].y()*spc.p[i].y() ); 
            std::vector<double> v = { radial };
            histChar.to_index(v);
            histCat[ v ] += 1;
            histChar[ v ] += atom[spc.p[i].id].charge;
          }
        }
        for ( auto i : agg ) {
          double radial = sqrt( spc.p[i].x()*spc.p[i].x() + spc.p[i].y()*spc.p[i].y() ); 
          std::vector<double> v = { radial };
          histChar.to_index(v);
          histChar[ v ] += atom[spc.p[i].id].charge;
        }
      }
    } // end of micro loop
  } // end of macro loop
   
  spc.save("state");               
  double scale = 1. / double(sample);
  histChar.save("histChar.dat", scale);
  histCat.save("histCat.dat", scale);
  cout << "total number of MC sweeps" << tot << " = " << sweeps << ", sampling events: " << sample << endl;
  cout << loop.info() + mv.info() + analysis.info() << endl;
}
