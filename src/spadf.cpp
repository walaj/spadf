#include <stdint.h>
#include <iostream>
#include <zlib.h>
#include <cassert>
#include <sstream>

#include <vector>
#include <cstring>

#include <memory>
#include <typeinfo>
#include <cmath>
#include <getopt.h>

#include <netcdf.h>
#include <Cell2.h>

///////
// MACROS
///////
#define NC_DEFINE_DIM(name, len, id) \
if (nc_def_dim(ncid, name, len, &id) != NC_NOERR) \
{ \
  std::cerr << "Error defining " << name << " dimension: " << nc_strerror(ncerr) << std::endl; \
  return 1; \
}

#define NC_DEFINE_VAR(name, type, rank, dimids, id) \
if (nc_def_var(ncid, name, type, rank, dimids, &id) != NC_NOERR) \
{ \
  std::cerr << "Error defining " << name << " variable: " << nc_strerror(ncerr) << std::endl; \
  return 1; \
}

#define NC_WRITE_VAR_FLOAT(name, id, data) \
if (nc_put_var_float(ncid, id, data) != NC_NOERR) \
{ \
  std::cerr << "Error writing " << name << " data: " << nc_strerror(ncerr) << std::endl; \
}

#define NC_WRITE_VAR_STRING(name, id, data) \
if (nc_put_var_string(ncid, id, data) != NC_NOERR) \
{ \
  std::cerr << "Error writing " << name << " data: " << nc_strerror(ncerr) << std::endl; \
}

#define FILE_NAME "test.nc"

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

  // file id
  int ncid;

  // dimension ids
  int markers_dimid, cells_dimid, coord_dimid, cell_meta_dimid; //, mnodes_dimid, snodes_dimid;
  const int markers_len = table.markers.size();
  const int meta_len = table.meta.size();  
  const int cells_len = table.num_cells();
  const int cell_meta_len = table.meta.size();
  const int coord_len = 2;

  // variable ids
  int marker_ids_varid, marker_names_varid, marker_meta_varid, cell_meta_varid;
  int cell_ids_varid, cells_areas_varid;
  int marker_table_varid, coord_table_varid, cell_meta_table_varid;
  //  int sgraph_varid, mgraph_varid;

  // create the netcdf file
  if (nc_create(FILE_NAME, NC_NETCDF4, &ncid) != NC_NOERR) {
    std::cerr << "Error creating netCDF file: " << nc_strerror(ncerr) << std::endl;
    return 1;
  }
  
  // define the dimensions
  NC_DEFINE_DIM("markers", markers_len, markers_dimid);
  NC_DEFINE_DIM("cells", cells_len, cells_dimid);
  NC_DEFINE_DIM("coordinates", coord_len, coord_dimid);
  NC_DEFINE_DIM("meta", cell_meta_len, cell_meta_dimid);  

  // define the variables
  int dimids_marker_table[] = {markers_dimid, cells_dimid};
  NC_DEFINE_VAR("marker_table", NC_FLOAT, 2, dimids_marker_table, marker_table_varid);
  int dimids_cell_meta_table[] = {cell_meta_dimid, cells_dimid};
  NC_DEFINE_VAR("meta_table", NC_FLOAT, 2, dimids_cell_meta_table, cell_meta_table_varid);
  
  NC_DEFINE_VAR("marker_names", NC_STRING, 1, &markers_dimid, marker_names_varid);
  NC_DEFINE_VAR("meta_names", NC_STRING, 1, &cell_meta_dimid, cell_meta_varid);  
    
  int dimids_coord_table[] = {cells_dimid, coord_dimid};
  NC_DEFINE_VAR("coord_table", NC_FLOAT, 2, dimids_coord_table, coord_table_varid);


  //////
  // add the data
  //////

  // Create a 2D array to hold the fluorescence data
  float* data_marker_table = (float*) malloc(markers_len * cells_len * sizeof(float));
  // Create a 2D array to hold the cell meta data  
  float* data_cellmeta_table = (float*) malloc(meta_len * cells_len * sizeof(float));  

  // Fill the 2D array with the data from the table object
  int i_meta = 0;
  int i_marker = 0;
  for (const auto& entry : table) {

    // Make sure the vector is the correct size
    if (entry.second.size() != cells_len) {
      std::cerr << "ERROR: the loaded quant table is ragged" << std::endl;
      return 1;
    }

    // don't add the x any y to either marker table or meta table
    if (entry.first == table.x || entry.first == table.y)
      continue;

    // Copy the marker values into the 2D array
    if (table.markers.count(entry.first)) {
      std::copy(entry.second.begin(), entry.second.end(), &data_marker_table[i_marker*cells_len]);
      i_marker++;
    }
    else if (table.meta.count(entry.first)) {
      std::copy(entry.second.begin(), entry.second.end(), &data_cellmeta_table[i_meta*cells_len]);
      i_meta++;
    } else {
      std::cerr << "ERROR. Element " << entry.first << " doesn't belong to markers or meta" << std::endl;
      return 1;
    }
  }    
    
    // Write the fluorescence table data to the netCDF file
  NC_WRITE_VAR_FLOAT("marker_table", marker_table_varid, data_marker_table);

  // free the marker table data
  delete data_marker_table;

  // Write the meta table data to the netCDF file
  NC_WRITE_VAR_FLOAT("meta_table", cell_meta_table_varid, data_cellmeta_table);
  
  // create the marker names data
  char** marker_array = new char*[markers_len];
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

  return 0;
  
}
