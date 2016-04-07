#include "MantidDataHandling/ImggAggregateWavelengths.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/make_unique.h"
#include "MantidKernel/MandatoryValidator.h"
#include "MantidKernel/StringTokenizer.h"

#include <fstream>

#include <Poco/File.h>
#include <Poco/DirectoryIterator.h>
#include <Poco/Path.h>

namespace Mantid {
namespace DataHandling {

const std::vector<std::string> ImggAggregateWavelengths::formatExtensionsShort{
    "fit"};
const std::vector<std::string> ImggAggregateWavelengths::formatExtensionsLong{
    "fits"};

const std::string ImggAggregateWavelengths::outPrefix = "bands_";
const std::string ImggAggregateWavelengths::indexRangesPrefix = "indices_";
const std::string ImggAggregateWavelengths::tofRangesPrefix = "tof_";

using Mantid::Kernel::Direction;
using Mantid::API::WorkspaceProperty;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ImggAggregateWavelengths)

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string ImggAggregateWavelengths::name() const {
  return "ImggAggregateWavelengths";
}

/// Algorithm's version for identification. @see Algorithm::version
int ImggAggregateWavelengths::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string ImggAggregateWavelengths::category() const {
  return "DataHandling\\Imaging";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string ImggAggregateWavelengths::summary() const {
  return "Aggregates images from multiple energy bands or wavelengths";
}

namespace {
const std::string PROP_INPUT_PATH = "InputPath";
const std::string PROP_OUTPUT_PATH = "OutputPath";
const std::string PROP_UNIFORM_BANDS = "UniformBands";
const std::string PROP_INDEX_RANGES = "IndexRanges";
const std::string PROP_TOF_RANGES = "ToFRanges";
const std::string PROP_OUTPUT_PREFIX = "OutputBandPrefix";
const std::string PROP_INPUT_FORMAT = "InputImageFormat";
const std::string PROP_OUTPUT_IMAGE_FORMAT = "OutputImageFormat";
const std::string PROP_OUTPUT_BIT_DEPTH = "OutputBitDepth";
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void ImggAggregateWavelengths::init() {
  declareProperty(Mantid::Kernel::make_unique<API::FileProperty>(
                      PROP_INPUT_PATH, "", API::FileProperty::Directory),
                  "The path (directory) to the input image files. It must "
                  "contain image files (radiography) or subdirectories "
                  "containing image files (radiography, one directory per "
                  "projection angle)");

  declareProperty(
      Kernel::make_unique<API::FileProperty>(PROP_OUTPUT_PATH, "",
                                             API::FileProperty::Directory),
      "The path (directory) where to generate the output image(s).");

  auto positive = boost::make_shared<Kernel::BoundedValidator<int>>();
  positive->setLower(0);
  declareProperty(
      PROP_UNIFORM_BANDS, 1, positive,
      "The number of output bands. The input bands are divided into uniform "
      "non-overlapping blocks and each of the output bands correspond to one "
      "block. This is a convenience particular case of the property " +
          PROP_INDEX_RANGES,
      Direction::Input);

  declareProperty(
      PROP_INDEX_RANGES, "",
      "A comma separated list of ranges of indices, where the "
      "indices refer to the input images, counting from 1. For "
      "example: 1-100, 101-200, 201-300. The number of ranges "
      "given will set the number of output bands. If you just need a "
      "certain"
      "number of bands split uniformly you can alternatively use "
      "the simple property " +
          PROP_UNIFORM_BANDS);

  // TODO: for later, when we have the required headers to do this
  declareProperty(
      PROP_TOF_RANGES, "",
      "A comma separated list of time-of-flight ranges given as for example "
      "5000-10000. These will be the boundaries in ToF that delimit the output "
      "bands. The units are as specified in the input images headers or "
      "metainformation (normally units of microseconds). The algorithm will "
      "produce as many output bands as ranges are "
      "given in this property. The ranges can overlap");

  declareProperty(
      PROP_OUTPUT_PREFIX, outPrefix,
      Kernel::make_unique<Kernel::MandatoryValidator<std::string>>(),
      "The output bands will be written into subdirectories with names where "
      "the prefix is prepended. In addition to this prefix, a second prefix "
      "that specifies whether the input bands were aggregated by indices or by "
      "time of flight is also appended. For example when running this "
      "algorithm using the property " +
          PROP_UNIFORM_BANDS + " or " + PROP_INDEX_RANGES +
          " the output subdirectoy names will look like '" + outPrefix +
          indexRangesPrefix + "1_1200' (where the 1 and 1200 derive from the "
                              "division into uniform non-overlaping blocks of "
                              "input bands, or the index ranges given. When "
                              "running the algorithm using the property " +
          PROP_TOF_RANGES + " sthe output subdirectory names will look like '" +
          outPrefix + tofRangesPrefix +
          "10000_50000' (where the 10000 and 50000 are the time of flight "
          "boundaries of the output band, using the same units as in the image "
          "headers).",
      Direction::Input);

  std::vector<std::string> imgFormat{"FITS"};
  declareProperty(
      PROP_INPUT_FORMAT, "FITS", "From the input directory(ies) use "
                                 "images in this format and ignore any "
                                 "other files",
      boost::make_shared<Mantid::Kernel::StringListValidator>(imgFormat),
      Direction::Input);

  declareProperty(
      PROP_OUTPUT_IMAGE_FORMAT, "FITS", "Produce output images in this format",
      boost::make_shared<Mantid::Kernel::StringListValidator>(imgFormat),
      Direction::Input);

  const std::vector<int> bitDepths{16};
  declareProperty(PROP_OUTPUT_BIT_DEPTH, 16,
                  boost::make_shared<Kernel::ListValidator<int>>(bitDepths),
                  "The bit depth or number of bits per pixel to use for the "
                  "output image(s). Only 16 bits is supported at the "
                  "moment.",
                  Direction::Input);
}

std::map<std::string, std::string> ImggAggregateWavelengths::validateInputs() {
  std::map<std::string, std::string> result;

  size_t optCount = 0;
  const int uniformBands = getProperty(PROP_UNIFORM_BANDS);
  if (uniformBands > 0)
    optCount++;
  const std::string indexRanges = getPropertyValue(PROP_INDEX_RANGES);
  if (!indexRanges.empty())
    optCount++;
  const std::string tofRanges = getPropertyValue(PROP_TOF_RANGES);
  if (!tofRanges.empty())
    optCount++;

  if (1 != optCount) {
    result[PROP_UNIFORM_BANDS] = "One and only one of the options " +
                                 PROP_UNIFORM_BANDS + ", " + PROP_INDEX_RANGES +
                                 ", " + PROP_TOF_RANGES + " has to be set.";
  }

  const std::string inputPath = getPropertyValue(PROP_INPUT_PATH);
  const std::string outputPath = getPropertyValue(PROP_OUTPUT_PATH);
  if (inputPath == outputPath)
    result[PROP_INPUT_PATH] = "The input and output paths should be different.";

  return result;
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void ImggAggregateWavelengths::exec() {
  const std::string inputPath = getPropertyValue(PROP_INPUT_PATH);
  const std::string outputPath = getPropertyValue(PROP_OUTPUT_PATH);
  const std::string indexRanges = getPropertyValue(PROP_INDEX_RANGES);
  const std::string tofRanges = getPropertyValue(PROP_TOF_RANGES);
  const int uniformBands = getProperty(PROP_UNIFORM_BANDS);

  if (uniformBands > 0)
    aggUniformBands(inputPath, outputPath, uniformBands);
  else if (!indexRanges.empty())
    aggIndexBands(inputPath, outputPath, indexRanges);
  else if (!tofRanges.empty())
    aggToFBands(inputPath, outputPath, tofRanges);
}

void ImggAggregateWavelengths::aggUniformBands(const std::string &inputPath,
                                               const std::string &outputPath,
                                               size_t bands) {

  Poco::Path path(inputPath);
  auto inputSubDirs = findInputSubdirs(path);

  if (inputSubDirs.empty()) {
    g_log.warning() << "Could not find any input files or directories in "
                    << inputPath
                    << ". Nothing will be produced in the output path."
                    << std::endl;
  }

  size_t count = 0;
  for (const auto &subdir : inputSubDirs) {
    processDirectory(subdir, bands, outputPath, outPrefix + indexRangesPrefix,
                     count++);
  }
}

void ImggAggregateWavelengths::aggIndexBands(const std::string &inputPath,
                                             const std::string &outputPath,
                                             const std::string &rangesSpec) {

  auto ranges = rangesFromStringProperty(rangesSpec, PROP_INDEX_RANGES);

  Poco::Path path(inputPath);
  auto inputSubDirs = findInputSubdirs(path);

  if (inputSubDirs.empty()) {
    g_log.warning() << "Could not find any input files or directories in "
                    << inputPath
                    << ". Nothing will be produced in the output path."
                    << std::endl;
  }

  size_t count = 0;
  for (const auto &subdir : inputSubDirs) {
    processDirectory(subdir, ranges, outputPath, outPrefix + indexRangesPrefix,
                     count++);
  }
}

void ImggAggregateWavelengths::aggToFBands(const std::string & /*inputPath*/,
                                           const std::string & /*outputPath*/,
                                           const std::string & /*ranges*/) {

  throw std::runtime_error("The property " + PROP_TOF_RANGES +
                           "is not supported at the moment");
}

std::vector<std::pair<size_t, size_t>>
ImggAggregateWavelengths::rangesFromStringProperty(
    const std::string &rangesSpec, const std::string &propName) {
  Mantid::Kernel::StringTokenizer rangeTokens(
      rangesSpec, ",", Mantid::Kernel::StringTokenizer::TOK_IGNORE_EMPTY |
                           Mantid::Kernel::StringTokenizer::TOK_TRIM);

  if (rangeTokens.count() <= 0) {
    throw std::invalid_argument(
        "Could not find any valid range of values in the property " + propName);
  }

  std::vector<std::pair<size_t, size_t>> ranges;
  for (const auto &str : rangeTokens) {
    const std::string sep = "-";
    Mantid::Kernel::StringTokenizer minMaxTokens(
        str, sep, Mantid::Kernel::StringTokenizer::TOK_IGNORE_EMPTY |
                      Mantid::Kernel::StringTokenizer::TOK_TRIM);
    if (2 != minMaxTokens.count()) {
      throw std::invalid_argument(
          "Could not parse a minimum and maximum value separated by '" + sep +
          "' from the string: " + str);
    }

    try {
      ranges.emplace_back(
          std::make_pair(boost::lexical_cast<size_t>(minMaxTokens[0]),
                         boost::lexical_cast<size_t>(minMaxTokens[1])));
    } catch (std::runtime_error &rexc) {
      throw std::runtime_error("Failed to parse this index range: " + str +
                               " . Details: " + std::string(rexc.what()));
    }
  }

  return ranges;
}

// just to compare two Poco::Path objects
struct PocoPathComp
    : public std::binary_function<Poco::Path, Poco::Path, bool> {
  bool operator()(const Poco::Path &lhs, const Poco::Path &rhs) const {
    return lhs.toString() < rhs.toString();
  }
};

bool ImggAggregateWavelengths::isSupportedExtension(
    const std::string &extShort, const std::string &extLong) {

  bool found = (formatExtensionsShort.end() !=
                std::find(formatExtensionsShort.cbegin(),
                          formatExtensionsShort.cend(), extShort)) ||
               (formatExtensionsLong.end() !=
                std::find(formatExtensionsLong.cbegin(),
                          formatExtensionsLong.cend(), extLong));
  return found;
}

/**
 * Find the (sub)directories that need processing. Each subdirectory
 * would normally correspond to a projection (angle). If an
 * (apparently) image file is found in the input path, it is assumed
 * that the input path is the only that needs processing (it should
 * contain images from a radiography or single projection) and
 * subdirectories are then ignored. Other files (for example .txt
 * files) are ignored in any case.
 *
 * @param path input path which may contain subdirectories
 *
 * @return subdirectories of the input path to process (or the input
 * path itself if it is the only one that will b processed).
 */
std::vector<Poco::Path>
ImggAggregateWavelengths::findInputSubdirs(const Poco::Path &path) {
  // Note: sorted directory iterators available only in poco >= 1.5.2
  std::vector<Poco::Path> dirsFound;

  Poco::File dir(path);
  if (!dir.isDirectory() || !dir.exists())
    return dirsFound;

  // could also use Poco::Glob to find the files as an alternative
  Poco::DirectoryIterator it(dir);
  Poco::DirectoryIterator end;
  while (it != end) {

    // there is at least one image file: take just the first level directory
    if (it->isFile()) {
      const std::string name = it.name();
      const std::string extShort = name.substr(name.size() - 3);
      const std::string extLong = name.substr(name.size() - 4);

      if (isSupportedExtension(extShort, extLong)) {
        // the input path contiains image files, take it
        dirsFound = {path};
        break;
      }
    }

    if (it->isDirectory()) {
      dirsFound.push_back(it.path());
    }

    ++it;
  }

  // this assumes the usual sorting of images of a stack (directory): a prefix,
  // and a sequence number (with a fixed number of digits).
  std::sort(dirsFound.begin(), dirsFound.end(), PocoPathComp());

  return dirsFound;
}

/**
 * Find images in a directory. It takes the (image) files that have
 * the expected/supported extension. Anything else (other files,
 * directories) is ignored.
 *
 * @param path directory where to look for images
 */
std::vector<Poco::Path>
ImggAggregateWavelengths::findInputImages(const Poco::Path &path) {
  // Note: sorted directory iterators available only in poco >= 1.5.2
  std::vector<Poco::Path> imgsFound;

  Poco::File dir(path);
  if (!dir.isDirectory() || !dir.exists())
    return imgsFound;

  Poco::DirectoryIterator it(dir);
  Poco::DirectoryIterator end;
  while (it != end) {

    if (it->isFile()) {
      const std::string summedSkipStr = "_SummedImg.";
      const std::string name = it.name();
      if (std::string::npos != name.find(summedSkipStr)) {
        ++it;
        continue;
      }

      const std::string expectedShort = "fit";
      const std::string expectedLong = "fits";
      const std::string extShort = name.substr(name.size() - 3);
      const std::string extLong = name.substr(name.size() - 4);
      if (isSupportedExtension(extShort, extLong)) {
        imgsFound.emplace_back(it.path());
      }
    }

    if (it->isDirectory()) {
      imgsFound.push_back(it.path());
    }

    ++it;
  }

  // this assumes the usual sorting of images of a stack (directory): a prefix,
  // and a sequence number (with a fixed number of digits).
  std::sort(imgsFound.begin(), imgsFound.end(), PocoPathComp());

  return imgsFound;
}

/**
 * Aggregate the images found in an input directory (can be all the
 * images from a radiography experiment, or all the images
 * corresponding to a single projection angle from a tomography
 * experiment).
 *
 * @param inDir where to find the input images.
 *
 * @param bands number of output bands, or number of bands into
 * which the input bands should be aggregated, non-overlapping and
 * uniformly distributed.
 *
 * @param outDir where to write the output image(s).
 *
 * @param prefix prefix for the image names. The indices will be
 * appended to the prefix
 */
void ImggAggregateWavelengths::processDirectory(const Poco::Path &inDir,
                                                size_t bands,
                                                const std::string &outDir,
                                                const std::string &prefix,
                                                size_t outImgIndex) {
  Mantid::API::MatrixWorkspace_sptr aggImg;

  size_t max = 0;
  // TODO: get max from directory (subdirs) listing
  auto images = findInputImages(inDir);

  if (images.empty()) {
    g_log.warning()
        << "Could not find any input image files in the subdirectory '"
        << inDir.toString() << "'. It will be ignored." << std::endl;
    return;
  }

  // TODO call the 'ranges' version with the convenient range
  const std::vector<std::pair<size_t, size_t>> ranges{std::make_pair(1, max)};
  processDirectory(inDir, ranges, outDir, prefix, outImgIndex);
}

API::MatrixWorkspace_sptr
ImggAggregateWavelengths::loadFITS(const Poco::Path &imgPath,
                                   const std::string &outName) {

  auto loader =
      Mantid::API::AlgorithmManager::Instance().createUnmanaged("LoadFITS");
  try {
    loader->initialize();
    loader->setPropertyValue("Filename", imgPath.toString());
    loader->setProperty("OutputWorkspace", outName);
    // this is way faster when loading into a MatrixWorkspace
    loader->setProperty("LoadAsRectImg", true);
  } catch (std::exception &e) {
    throw std::runtime_error("Failed to initialize the algorithm to "
                             "load images. Error description: " +
                             std::string(e.what()));
  }

  try {
    loader->execute();
  } catch (std::exception &e) {
    throw std::runtime_error(
        "Failed to load image. Could not load this file as a "
        "FITS image: " +
        std::string(e.what()));
  }

  if (!loader->isExecuted()) {
    throw std::runtime_error(
        "Failed to load image correctly. Note that even though "
        "the image file has been loaded it seems to contain errors.");
  }

  Mantid::API::MatrixWorkspace_sptr ws;
  try {
    Mantid::API::WorkspaceGroup_sptr wsg;
    const auto &ads = Mantid::API::AnalysisDataService::Instance();
    wsg = ads.retrieveWS<Mantid::API::WorkspaceGroup>(outName);
    ws = ads.retrieveWS<Mantid::API::MatrixWorkspace>(
        wsg->getNames()[wsg->size() - 1]);
  } catch (std::exception &e) {
    throw std::runtime_error(
        "Could not load image contents for file '" + imgPath.toString() +
        "'. An unrecoverable error "
        "happened when trying to load the image contents. Cannot "
        "display it. Error details: " +
        std::string(e.what()));
  }

  return ws;
}

/**
 * Add a workspace into an accumulator workspace ('toAdd' into
 * 'accum'). This method blatantly ignores the X values and the E
 * (error) values. We're only interested in the Y values (pixels,
 * counts).
 *
 * @param accum workspace where to add a new workspace
 * @param toAdd workspace to add
 */
void ImggAggregateWavelengths::aggImage(API::MatrixWorkspace_sptr accum,
                                        const API::MatrixWorkspace_sptr toAdd) {
  const size_t sizeX = accum->blocksize();
  const size_t sizeY = accum->getNumberHistograms();

  const size_t sizeXIn = accum->blocksize();
  const size_t sizeYIn = accum->getNumberHistograms();

  if (sizeX != sizeXIn || sizeY != sizeYIn) {
    throw std::runtime_error(
        "Cannot add images of different dimensions. The first image has "
        "dimensions " +
        std::to_string(sizeY) + " rows by " + std::to_string(sizeX) +
        " columns whereas the second image has dimensions " +
        std::to_string(sizeYIn) + " rows by " + std::to_string(sizeXIn) +
        " columns.");
  }

  for (size_t row = 0; row < sizeY; row++) {
    Mantid::API::ISpectrum *spectrum = accum->getSpectrum(row);
    Mantid::API::ISpectrum *specIn = toAdd->getSpectrum(row);
    auto &dataY = spectrum->dataY();
    const auto &dataYIn = specIn->readY();
    std::transform(dataY.begin(), dataY.end(), dataYIn.cbegin(), dataY.begin(),
                   std::plus<double>());
    // for (size_t col = 0; col < sizeX; col++) {
    //  dataY[col] += dataYIn[col];
    //}
  }
}

void ImggAggregateWavelengths::saveAggImage(
    const API::MatrixWorkspace_sptr accum, const std::string &outDir,
    const std::string &prefix, size_t outImgIndex) {
  // for example 'bands_sum_00030'
  std::ostringstream sstr;
  sstr << std::setw(5) << std::setfill('0') << outImgIndex;
  const std::string outName = prefix + "sum_" + sstr.str();

  Poco::Path outPath(outDir);
  Poco::File dirFile(outPath);
  if (!dirFile.isDirectory() || !dirFile.exists()) {
    g_log.information() << "Cannot save output image into '"
                        << outPath.toString()
                        << "'. It is not an existing directory." << std::endl;
    return;
  }

  outPath.append(outName);
  // only FITS support for now
  const std::string fullName = outPath.toString();
  saveFITS(accum, fullName);
  g_log.information() << "Saved output aggregated image into: " << fullName
                      << std::endl;
}

namespace {
// minimalistic. At a very least we should add ToF, time bin, counts, triggers,
// etc.
const size_t maxLenHdr = 80;
const std::string FITSHdrFirst =
    "SIMPLE  =                    T / file does conform to FITS standard";
const std::string FITSHdrBitDepth =
    "BITPIX  =                   16 / number of bits per data pixel";
const std::string FITSHdrAxes =
    "NAXIS   =                    2 / number of data axes";
const std::string FITSHdrExtensions =
    "EXTEND  =                    T / FITS dataset may contain extensions";
const std::string FITSHdrRefComment1 = "COMMENT   FITS (Flexible Image "
                                       "Transport System) format is defined in "
                                       "'Astronomy";
const std::string FITSHdrRefComment2 =
    "COMMENT   and Astrophysics', volume 376, page 359; bibcode: "
    "2001A&A...376..359H";
const std::string FITSHdrEnd = "END";

void writeFITSHeaderEntry(const std::string &hdr, std::ofstream &file) {
  static const std::array<char, maxLenHdr> zeros{};

  size_t count = hdr.size();
  if (count >= maxLenHdr)
    count = maxLenHdr;

  file.write(hdr.c_str(), sizeof(char) * count);
  file.write(&zeros[0], maxLenHdr - count);
}

void writeFITSHeaderAxesSizes(const API::MatrixWorkspace_sptr img,
                              std::ofstream &file) {
  const std::string sizeX = std::to_string(img->blocksize());
  const std::string sizeY = std::to_string(img->getNumberHistograms());

  const size_t fieldWidth = 20;
  std::stringstream axis1;
  axis1 << "NAXIS1  = " << std::setw(fieldWidth) << sizeX
        << " / length of data axis 1";
  writeFITSHeaderEntry(axis1.str(), file);

  std::stringstream axis2;
  axis2 << "NAXIS2  = " << std::setw(fieldWidth) << sizeY
        << " / length of data axis 2";
  writeFITSHeaderEntry(axis2.str(), file);
}

// FITS headers consist of subblocks of 36 entries/lines. This method
// is to write the "padding" lines required to have 36 in a block
void writePaddingFITSHeaders(size_t count, std::ofstream &file) {
  static const std::array<char, maxLenHdr> zeros{};

  for (size_t i = 0; i < count; ++i) {
    file.write(&zeros[0], maxLenHdr);
  }
}

void writeFITSHeaderBlock(const API::MatrixWorkspace_sptr img,
                          std::ofstream &file) {
  writeFITSHeaderEntry(FITSHdrFirst, file);
  writeFITSHeaderEntry(FITSHdrBitDepth, file);
  writeFITSHeaderEntry(FITSHdrAxes, file);
  writeFITSHeaderAxesSizes(img, file);
  writeFITSHeaderEntry(FITSHdrExtensions, file);
  writeFITSHeaderEntry(FITSHdrRefComment1, file);
  writeFITSHeaderEntry(FITSHdrRefComment2, file);
  writeFITSHeaderEntry(FITSHdrEnd, file);

  const size_t entriesPerHDU = 36;
  writePaddingFITSHeaders(entriesPerHDU - 9, file);
}

void writeFITSImageMatrix(const API::MatrixWorkspace_sptr img,
                          std::ofstream &file) {
  const size_t sizeX = img->blocksize();
  const size_t sizeY = img->getNumberHistograms();

  for (size_t row = 0; row < sizeY; row++) {
    Mantid::API::ISpectrum *spectrum = img->getSpectrum(row);
    const auto &dataY = spectrum->readY();
    for (size_t col = 0; col < sizeX; col++) {
      const size_t bytespp = 2;
      int16_t pixelVal = static_cast<uint16_t>(dataY[col]);

      // change endianness: to sequence of bytes in big-endian
      uint8_t bytesPixel[bytespp];
      uint8_t *iter = reinterpret_cast<uint8_t *>(&pixelVal);
      std::reverse_copy(iter, iter + bytespp, bytesPixel);

      file.write(reinterpret_cast<const char *>(&bytesPixel[0]),
                 sizeof(bytesPixel));
    }
  }
}
}

// TODO: this is a very early version of what should become an
// algorithm of its own
void ImggAggregateWavelengths::saveFITS(const API::MatrixWorkspace_sptr accum,
                                        const std::string &filename) {
  std::ofstream outfile(filename + ".fits", std::ofstream::binary);

  writeFITSHeaderBlock(accum, outfile);
  writeFITSImageMatrix(accum, outfile);
}

void ImggAggregateWavelengths::processDirectory(
    const Poco::Path &inDir,
    const std::vector<std::pair<size_t, size_t>> &ranges,
    const std::string outDir, const std::string prefix, size_t outImgIndex) {

  auto imgFiles = findInputImages(inDir);

  auto it = std::begin(imgFiles);
  const std::string wsName = "__ImggAggregateWavelengths_fits_seq";
  const std::string wsNameFirst = wsName + "_first";
  API::MatrixWorkspace_sptr accum = loadFITS(*it, wsNameFirst);
  ++it;

  for (auto end = std::end(imgFiles); it != end; ++it) {
    // load one more
    const API::MatrixWorkspace_sptr img = loadFITS(*it, wsName);

    // add
    aggImage(accum, img);

    // clear image/workspace. TODO: This is a big waste of
    // allocations/deallocations
    Mantid::API::AnalysisDataService::Instance().remove(wsName);
  }
  // save
  saveAggImage(accum, outDir, prefix, outImgIndex);
  Mantid::API::AnalysisDataService::Instance().remove(wsNameFirst);
}

} // namespace DataHandling
} // namespace Mantid
