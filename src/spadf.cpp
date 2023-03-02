#include <stdint.h>
#include <iostream>
#include <zlib.h>
#include <cassert>
#include <sstream>

#include <vector>
#include <cstring>

#include "spadf.h"
#include "spadf_utils.h"

#include <memory>
#include <typeinfo>
#include <cmath>
#include <getopt.h>

#include <netcdf.h>
#include <Cell2.h>

#include "umappp/Umap.hpp"
#include "umappp/NeighborList.hpp"
#include "umappp/combine_neighbor_sets.hpp"
#include "umappp/find_ab.hpp"
#include "umappp/neighbor_similarities.hpp"
#include "umappp/optimize_layout.hpp"
#include "umappp/spectral_init.hpp"
#include "knncolle/knncolle.hpp"

namespace opt {
  static bool verbose = false;
  static std::string infile;
  static std::string quantfile;
  static std::string markerfile;  
  static std::string outfile;
  static std::string module;

  static std::string redfile;
  static std::string greenfile;
  static std::string bluefile;
  static int threads = 1;
}

#define DEBUG(x) std::cerr << #x << " = " << (x) << std::endl

#define TVERB(msg) \
  if (opt::verbose) {		   \
    std::cerr << msg << std::endl; \
  }

static const char* shortopts = "hvr:g:b:q:c:m:";
static const struct option longopts[] = {
  { "verbose",                    no_argument, NULL, 'v' },
  { "marker-file",                required_argument, NULL, 'm' },  
  { "threads",                    required_argument, NULL, 'c' },
  { "quant-file",                 required_argument, NULL, 'q' },  
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: spadf [module] <options> \n"
"Modules:\n"
"  debug - sandbox\n"
  "\n";
  
static int debugfunc();
static void parseRunOptions(int argc, char** argv);

int main(int argc, char **argv) {
  
  // Check if a command line argument was provided
  if (argc < 2) {
    std::cerr << "Error: missing command line argument" << std::endl;
    std::cerr << RUN_USAGE_MESSAGE;
    return 1;
  }

  parseRunOptions(argc, argv);
  
  // get the module
  if (opt::module == "debug") {
    return(debugfunc());
  } else {
    assert(false);
  }
  
  return 1;
}

// parse the command line options
static void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 1) 
    die = true;

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'c' : arg >> opt::threads; break;
    case 'm' : arg >> opt::markerfile; break;            
    case 'q' : arg >> opt::quantfile; break;      
    case 'h' : help = true; break;
    default: die = true;
    }
  }
  
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::module.empty()) {
      opt::module = argv[optind];      
    } 
    optind++;
  }

  if (! (opt::module == "debug")) {
    std::cerr << "Module " << opt::module << " not implemented" << std::endl;
    die = true;
  }
  
  if (die || help) 
    {
      std::cerr << "\n" << RUN_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

static int debugfunc() {

  ////////
  // create
  ////////
  CellTable table(opt::quantfile.c_str(), opt::markerfile.c_str());

  spadf spa("test2.nc", table);

  spa.subsample_cells(2, 1337);
  
  spa.knn_marker();
  
  spa.Close();

  return 1;
  
  // create the marker names data
  /*  char** marker_array = new char*[markers_len];
  table.FillMarkerNames(marker_array);

  // write marker strings
  const char** marker_array_const = const_cast<const char**>(marker_array);
  NC_WRITE_VAR_STRING("marker_names", marker_names_varid, marker_array_const);
  
  // create the meta names data
  char** meta_array = new char*[meta_len];
  table.FillMetaNames(meta_array);  
  
  // write some strings
  const char** meta_array_const = const_cast<const char**>(meta_array);
  NC_WRITE_VAR_STRING("meta_names", cell_meta_varid, meta_array_const);
  
  // Free the memory allocated for each marker name
  for (size_t i = 0; i < markers_len; i++)
    delete[] marker_array[i];
  delete[] marker_array;
  
  for (size_t i = 0; i < meta_len; i++)
    delete[] meta_array[i];
  delete[] meta_array;
  
  // Create a 2D array to hold the coordinate data
  float* data_coord_table = (float*) malloc(coord_len * cells_len * sizeof(float));
  
  // Fill the 2D array with the data from the table object
  std::copy(table.GetCoord('x').begin(), table.GetCoord('x').end(), data_coord_table);
  std::copy(table.GetCoord('y').begin(), table.GetCoord('y').end(), &data_coord_table[cells_len]);  

  // Write the coord table data to the netCDF file
  NC_WRITE_VAR_FLOAT("coord_table", coord_table_varid, data_coord_table);

  delete data_coord_table;
  
  // close the netcdf file
  if (nc_close(ncid) != NC_NOERR) {
    std::cerr << "Error closing netCDF file: " << nc_strerror(ncerr) << std::endl;
    return 1;
  }
  */
  
  return 0;
  
}

int spadf::get_marker_rows(size_t row, size_t nrows, float* data) {

  // NEED error checking etc
  
  // have to specify the corner to start and the number of rows
  size_t startp[2] = {0,row};
  size_t countp[2] = {markers_len, nrows};  

  int status = nc_get_vara_float(ncid,
				 marker_table_varid,
				 startp,
				 countp,
				 data);

  return status;
  
}

spadf::spadf(const char* file, const CellTable& table) {

  // create the netcdf file
  if (nc_create(file, NC_NETCDF4, &ncid) != NC_NOERR) {
    std::cerr << "Error creating netCDF file: " << nc_strerror(ncerr) << std::endl;
    return;
  }

  // define the dimensions
  markers_len = table.markers.size();

  if (!markers_len) {
    std::cerr << "Warning: Markers len is 0, no data to be stored." << std::endl;
  }

  cells_len = table.num_cells();
  if (!cells_len) {
    std::cerr << "Warning: Cells len is 0, no data to be stored." << std::endl;
  }

  meta_cells_len = table.meta.size();  
  coords_len = 2;

  // define the dimensions
  NC_DEFINE_DIM("markers", markers_len, markers_dimid);
  NC_DEFINE_DIM("cells", cells_len, cells_dimid);
  NC_DEFINE_DIM("coordinates", coords_len, coords_dimid);
  NC_DEFINE_DIM("meta_markers", meta_markers_len, meta_markers_dimid);
  NC_DEFINE_DIM("meta_cells", meta_cells_len, meta_cells_dimid);    

  ////////////////
  // define the variables
  ////////////////

  //// 1D arrays
  
  // marker array of marker names
  NC_DEFINE_VAR("marker_names", NC_STRING, 1, &markers_dimid, markers_varid);
  
  // cell array of cell names
  NC_DEFINE_VAR("cell_names", NC_STRING, 1, &cells_dimid, cells_varid);  

  // cella meta array
  NC_DEFINE_VAR("meta_cells", NC_STRING, 1, &meta_cells_dimid, meta_cells_varid);
  
  //// 2D tables

  // marker-cell table of fluorescence values
  int dimids_marker_table[] = {markers_dimid, cells_dimid};
  NC_DEFINE_VAR("marker_table", NC_FLOAT, 2, dimids_marker_table, marker_table_varid);

  // cell-coordinate table
  int dimids_coord_table[] = {cells_dimid, coords_dimid};
  NC_DEFINE_VAR("coord_table", NC_FLOAT, 2, dimids_coord_table, coord_table_varid);
  
  // cell-meta table
  int dimids_meta_table[] = {cells_dimid, meta_cells_dimid};
  NC_DEFINE_VAR("meta_cells_table", NC_FLOAT, 2, dimids_meta_table, cells_meta_table_varid);

  ////////
  // add the data
  ///////

  if (AddMarkerTable(table) < 0)
    return;

  if (AddMetaTable(table) < 0)
    return;

  // write marker strings
  char** marker_array = new char*[markers_len];
  table.FillMarkerNames(marker_array);
  const char** marker_array_const = const_cast<const char**>(marker_array);
  NC_WRITE_VAR_STRING("marker_names", markers_varid, marker_array_const);
  for (size_t i = 0; i < markers_len; i++)
    delete[] marker_array[i];
  delete[] marker_array;

  // write meta strings
  char** meta_array = new char*[meta_cells_len];
  table.FillMarkerNames(meta_array);
  const char** meta_array_const = const_cast<const char**>(meta_array);
  NC_WRITE_VAR_STRING("meta_cells", meta_cells_varid, meta_array_const);
  for (size_t i = 0; i < meta_cells_len; i++)
    delete[] meta_array[i];
  delete[] meta_array;

}

int spadf::subsample_cells(int n, int seed) {

  // make sure we are subsampling less than original
  if (n > cells_len) {
    std::cerr << "Error: trying to sample " << n << " cells from " << cells_len << std::endl;
    return -1;
  }

  // create an inclusion set
  std::set<int> inc_set = sampleWithoutReplacement(n, cells_len, seed);

  // allocate the new table;
  float* new_marker_table = (float*) malloc(markers_len * n * sizeof(float)); 

  // temp holder for data row
  float* f = (float*) malloc(markers_len * sizeof(float));

  // subsample the cells
  int k = 0; // track how many 
  for (int i = 0; i < cells_len; i++) {

    // if not included in the random inclusion set, skip
    if (!inc_set.count(i))
      continue;

    std::cerr << " i " << i << std::endl;

    // get a single row at position (cell) i
    get_marker_rows(i, 1, f);

    // copy it over 
    memcpy(new_marker_table + k*markers_len, f, markers_len*sizeof(float));
    k++;
    
  }

  return 0;
  
}

int spadf::knn_marker() {

  // 2 comes from two output axes of umap
  std::vector<double> embedding(cells_len * 2);

  umappp::Umap x;
  x.set_num_threads(opt::threads);
  //x.set_num_neighbors();

  float* data= (float*) malloc(cells_len * markers_len * sizeof(float));
  nc_get_var_float(ncid, marker_table_varid, data);

  for (int i = 0; i < 4; i++)
    std::cerr << data[i] << std::endl;
  
  //std::vector<double> data = {}; //table.ColumnMajor();

  int ndim = markers_len;
  int nobs = cells_len;
    
  // find K nearest neighbors
  std::cerr << "...finding K nearest-neighbors" << std::endl;  
  int num_neighbors = 300;

  knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, data);
  std::cerr << " done initializeing " << std::endl;
  const size_t N = nobs; //searcher->nobs();
  umappp::NeighborList<float> output(N);
  std::cerr << " opt threads " << opt::threads << std::endl;
  #pragma omp parallel for num_threads(opt::threads) //rparams.nthreads)
  for (size_t i = 0; i < N; ++i) {
    if (i % 20000 == 0)
      std::cerr << " i " << AddCommas(i) << " thread " << omp_get_thread_num() << " K " << num_neighbors << std::endl;
    output[i] = searcher.find_nearest_neighbors(i, num_neighbors);
  }
  
  return 0;
  
}

int spadf::_table_check(const CellTable& table) {

  // check that the input table is filled in
  if (table.table.size() == 0) {
    std::cerr << "Warning: Cell table is empty." << std::endl;
    return -1;
  }

  // check that the input table is not ragged
  int data_len = -1;
  for (const auto& e : table.table) {
    // store the length of the first array
    if (data_len == -1) {
      data_len = e.second.size();
      continue;
    }

    if (e.second.size() != data_len) {
      std::cerr << "Warning: ragged array of cell table" << std::endl;
      return -1;
    }
  }

  // check that the data length is expected
  if (data_len != cells_len) {
    std::cerr << "Warning: Number of cells in table " << data_len <<
      " is different than expected " << cells_len << std::endl;
    return -1;
  }

  // check that the expected marker length is the same
  if (table.markers.size() != markers_len) {
    std::cerr << "Warning: Number of markers in table " << table.markers.size() <<
      " is different than expected " << markers_len << std::endl;
    return -1;
  }

  // check that the expected cell meta length is the same
  if (table.meta.size() != meta_cells_len) {
    std::cerr << "Warning: Number of cell meta columns in table " << table.meta.size() <<
      " is different than expected " << meta_cells_len << std::endl;
    return -1;
  }

  // check that the cell table is not ragged
  for (const auto& entry : table.table) {
    
    // Make sure the vector is the correct size
    if (entry.second.size() != cells_len) {
      std::cerr << "ERROR: the loaded quant table is ragged" << std::endl;
      return -1;
    }
  }
    
  return data_len;
  
}

int spadf::Close() {

  // close the netcdf file
  if (nc_close(ncid) != NC_NOERR) {
    std::cerr << "Error closing netCDF file: " << nc_strerror(ncerr) << std::endl;
    return 1;
  }

  return 0;
}

int spadf::AddMarkerTable(const CellTable& table) {

  int data_len = _table_check(table);
  
  // Create a 2D array to hold the fluorescence data
  //std::vector<float> *data_table = new std::vector<float>(markers_len * cells_len);
  float* data_table =
    (float*)malloc(markers_len * cells_len * sizeof(float));
  
  //std::unique_ptr<float*> data_table =
  //  reinterpret_cast<float*>(malloc(markers_len * cells_len * sizeof(float)));

  int i = 0;
  for (const auto& entry : table) {
    
    // don't add the x any y to either marker table or meta table
    if (entry.first == table.x || entry.first == table.y)
      continue;
    
    // Copy the marker values into the 2D array
    if (table.markers.count(entry.first)) {
      std::copy(entry.second.begin(), entry.second.end(), &data_table[i*cells_len]);
      i++;
    }
  }

  // check that the correct number of markers were added
  if (i != markers_len) {
    std::cerr << "Warning: Added " << i <<
      " markers but expected " << markers_len << std::endl;
    free(data_table);
    return -1;
  }

  // Write the fluorescence table data to the netCDF file
  NC_WRITE_VAR_FLOAT("marker_table", marker_table_varid, data_table);

  free(data_table);
  return 0;
}

int spadf::AddMetaTable(const CellTable& table) {

  int data_len = _table_check(table);
  if (data_len < 0) // error already thrown
    return -1;
  
  // Create a 2D array to hold the fluorescence data
  float* data_table =
    (float*)malloc(markers_len * cells_len * sizeof(float));
  
  int i = 0;
  for (const auto& entry : table.table) {
    
    // don't add the x any y to either marker table or meta table
    if (entry.first == table.x || entry.first == table.y)
      continue;
    
    // Copy the marker values into the 2D array
    if (table.meta.count(entry.first)) {
      std::copy(entry.second.begin(), entry.second.end(), &data_table[i*cells_len]);
      i++;
    }
  }

  // Check that the correct number of markers were added
  if (i != meta_cells_len) {
    std::cerr << "Warning: Added " << i <<
      " markers but expected " << meta_cells_len << std::endl;
    free(data_table);
    return -1;
  }

  // Write the fluorescence table data to the netCDF file
  NC_WRITE_VAR_FLOAT("meta_cells_table", cells_meta_table_varid, data_table);

  free(data_table);
  return 0;
  
}
