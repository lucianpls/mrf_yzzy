#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <gdal.h>
#include <cpl_string.h>

using namespace std;

int Usage(const char *message = nullptr, int retcode = 1) {
    if (message)
        cerr << message << endl;

    cerr << "mrf_yzzy transposes the data in a 3rD MRF by swapping the Y and Z axis" << endl
        << "Usage:" << endl
        << "mrf_yzzy [-z ZPageSize] [-v] [-g] in.mrf out.mrf" << endl << endl
        << "\t-z ZPageSize : Set the output Y pagesize" << endl
        << "\t-v : verbose" << endl
        << "\t-g copies the input projection and the area info, which will be wrong anyhow" << endl;

    return retcode;
}

int main(int argc, char **argv) {
    bool verbose = false;
    // Preserve the input geoprojection
    bool geo = false;
    int psz = 0; // No default
    GDALAllRegister();

    GDALDriverH d_mrf = GDALGetDriverByName("MRF");
    if (!d_mrf)
        return Usage("MRF driver not found");

    std::vector<std::string> fnames;

    // Pick up the GDAL options
    int nArgc = GDALGeneralCmdLineProcessor(argc, &argv, 0);
    if (nArgc < 1)
        exit(-nArgc);

    for (int iArg = 1; iArg < nArgc; iArg++)
    {
        if (EQUAL(argv[iArg], "-z")) {
            psz = atoi(argv[++iArg]);
        }
        else if (EQUAL(argv[iArg], "-v")) {
            verbose = true;
        }
        else if (EQUAL(argv[iArg], "-g")) {
            geo = true;
        } else
            fnames.push_back(argv[iArg]);
    }

    if (fnames.size() != 2)
        return Usage();

    string SourceName(fnames[0]), TargetName(fnames[1]);

    CPLPushErrorHandler(CPLQuietErrorHandler);
    GDALDatasetH hDatasetin = GDALOpen(SourceName.c_str(), GA_ReadOnly);
    CPLPopErrorHandler();

    if (hDatasetin == NULL) {
        CPLError(CE_Failure, CPLE_AppDefined, "Can't open source file %s for reading", SourceName.c_str());
        return 0;
    }

    if (!EQUAL(GDALGetDriverShortName(GDALGetDatasetDriver(hDatasetin)), "MRF"))
        return Usage("Input file is not MRF", 2);

    double gt[6] = {0, 0, 0, 0, 0, 0};
    CPLString projection(GDALGetProjectionRef(hDatasetin));
    char **md = GDALGetMetadata(hDatasetin, "IMAGE_STRUCTURE");
    if (!CSLFetchNameValue(md, "ZSIZE"))
        return Usage("Source is not a 3-rd dimension MRF", 2);
    int zsz = atoi(CSLFetchNameValue(md, "ZSIZE"));
    int csz = GDALGetRasterCount(hDatasetin);
    GDALRasterBandH b1 = GDALGetRasterBand(hDatasetin, 1);
    int xsz = GDALGetRasterBandXSize(b1);
    int ysz = GDALGetRasterBandYSize(b1);

    // Get the source geotransform and convert it for the output, preserving the area
    GDALGetGeoTransform(hDatasetin, gt);
    // gt[5] is the new y resolution, should be adjusted based on the new Y dimension
    gt[5] *= double(ysz) / double(zsz);

    // The NoData and Min-Max should be done per band instead of the current solution
    int bHasNoData = false;
    double nd = GDALGetRasterNoDataValue(b1, &bHasNoData);
    int bHasStats = false;

    // Get Stats if present
    double min_v, max_v, mean_v, stdd_v;
    bHasStats = (CE_None == GDALGetRasterStatistics(b1, TRUE, FALSE, &min_v, &max_v, &mean_v, &stdd_v));

    int pszx, pszy;
    GDALGetBlockSize(b1, &pszx, &pszy);

    GDALDataType dt = GDALGetRasterDataType(b1);
    int dtsz = GDALGetDataTypeSizeBytes(dt);

    // Checks and adjustments
    if (!psz)
        psz = pszx;

    char **copt = NULL;
    char **freeopt = NULL;

    while (md && *md) {
//        cout << *md << endl;
        if (STARTS_WITH_CI(*md, "COMPRESSION=")) {
            copt = CSLAppendPrintf(copt, "COMPRESS=%s", strstr(*md, "=") + 1);;
        }
        else if (STARTS_WITH_CI(*md, "ZSLICE=")
            || STARTS_WITH_CI(*md, "ZSIZE=")
            || STARTS_WITH_CI(*md, "V2=")
            )
        {
            // Removed, modified or ignored
        }
        // Free options
        else if (STARTS_WITH_CI(*md, "V1=")
            || STARTS_WITH_CI(*md, "GZ=")
            || STARTS_WITH_CI(*md, "ZSTD=")
            || STARTS_WITH_CI(*md, "RAWZ=")
            || STARTS_WITH_CI(*md, "DEFLATE=")
            || STARTS_WITH_CI(*md, "LERC_PREC=")
            )
        {
            freeopt = CSLAddString(freeopt, *md);
        }
        else { // Have no idea, assume create option
            copt = CSLAddString(copt, *md);
        }
        md++;
    }

    // Set the free form options, if any
    CPLString fopt;
    for (int i = 0; i < CSLCount(freeopt); i++) {
        if (fopt.size())
            fopt += " ";
        fopt += freeopt[i];
    }

    if (fopt.size())
        copt = CSLAddNameValue(copt, "OPTIONS", fopt);
    CSLDestroy(freeopt);

    // Add the known options
    copt = CSLAppendPrintf(copt, "BLOCKXSIZE=%d", pszx);
    copt = CSLAppendPrintf(copt, "BLOCKYSIZE=%d", psz);
    copt = CSLAppendPrintf(copt, "ZSIZE=%d", ysz);

    if (verbose) {
        md = copt;
        while (md && *md) {
            cout << *md << endl;
            md++;
        }
    }

    GDALClose(hDatasetin);

    // Operating on a block of size
    size_t BSZ = static_cast<size_t>(csz) * psz * pszy * pszx * dtsz;

    // These are the input strides
    int pix_stride = dtsz;
    int line_stride = pszx * pix_stride;
    int z_stride = pszy * line_stride;
    int band_stride = psz * z_stride;

    //char *outbuffer = reinterpret_cast<char *>(malloc(BSZ));
    char *buffer  = reinterpret_cast<char *>(malloc(BSZ));

    if (!buffer)
        return Usage(CPLOPrintf("Failed to allocate buffer of size %llu", BSZ), 3);
    if (verbose)
        cout << "Using an " << BSZ << " sized buffer\n";

    // Reading, Loop over y, x, z and c.
    // Start refers to input
    // End refers to output

    for (int startz = 0; startz < zsz; startz += psz) {
        int dz = min(psz, zsz - startz);

        vector<GDALDatasetH> inh(dz);
        for (int z = 0; z < dz; z++) {
            CPLString SName;
            SName.Printf("%s:MRF:Z%d", SourceName.c_str(), startz + z);
            GDALDatasetH  h = GDALOpen(SName, GA_ReadOnly);
            inh[z] = h;
        }

        for (int starty = 0; starty < ysz; starty += pszy) {
            int dy = min(pszy, ysz - starty);

            vector<GDALDatasetH> outh(dy);
            for (int z = 0; z < dy; z++) {
                CPLString DName;
                DName.Printf("%s:MRF:Z%d", TargetName.c_str(), starty + z);
                GDALDatasetH h = GDALCreate(d_mrf, DName.c_str(), xsz, zsz, csz, dt, copt);
                GDALRasterBandH b = GDALGetRasterBand(h, 1);
                if (bHasNoData)
                    GDALSetRasterNoDataValue(b, nd);
                if (bHasStats)
                    GDALSetRasterStatistics(b, min_v, max_v, mean_v, stdd_v);
                if (geo) {
                    GDALSetProjection(h, projection);
                    GDALSetGeoTransform(h, gt);
                }
                outh[z] = h;
            }

            // This loop does a full cublock
            for (int startx = 0; startx < xsz; startx += pszx) {
                int dx = min(pszx, xsz - starty);
                cout << "Processing " << startx << "," << starty << "," << startz << endl;
                // fprintf(stderr, "Processing %d,%d,%d\n", startx, starty, startz);

                // Read a cublock
                for (int z = 0; z < dz; z++) {
                    // Each one is considered a different dataset
                    hDatasetin = inh[z];
                    //fprintf(stderr,
                    //    "Reading Z%d %d,%d - %d,%d %d stride %d %d %d\n",
                    //    startz + z, startx, starty, dx, dy,
                    //    z_stride * (z + startz), pix_stride, line_stride, band_stride
                    //);
                    GDALDatasetRasterIO(hDatasetin, GF_Read,
                        startx, starty, dx, dy,
                        buffer + z_stride * z, dx, dy,
                        dt, csz, NULL,
                        pix_stride, line_stride, band_stride
                    );
                }

                // Write a cublock
                for (int endz = 0; endz < dy; endz++) {
                    GDALDatasetH hDatasetout = outh[endz];
                    //fprintf(stderr,
                    //    "Writing Z%d %d,%d - %d,%d %d stride %d %d %d\n",
                    //    starty + endz, startx, startz, dx, dy,
                    //    endz * line_stride, pix_stride, z_stride, band_stride
                    //);
                    GDALDatasetRasterIO(hDatasetout, GF_Write,
                        startx, startz, dx, dz,
                        buffer + endz * line_stride, dx, dz,
                        dt, csz, NULL,
                        pix_stride, z_stride, band_stride
                    );
                }
            }

            for (int i = 0; i < outh.size(); i++)
                GDALClose(outh[i]);
        }
        for (int i = 0; i < inh.size(); i++)
            GDALClose(inh[i]);
    }

    CSLDestroy(copt);
    return 0;
}