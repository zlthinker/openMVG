//
// Created by Larry CHOU on 15/12/12.
//

#include "third_party/easyexif/exif.h"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "nonFree/sift/SIFT_describer.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>
#include <dirent.h>
#include <third_party/ceres-solver/internal/ceres/gtest/gtest.h>
#include <queue>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;
using namespace svg;
using namespace std;

/// Read intrinsic K matrix from a file (ASCII)
/// F 0 ppx
/// 0 F ppy
/// 0 0 1
int readIntrinsic(const std::string & fileName, Mat3 & K);
double zncc(Image<unsigned char> image1, int x1, int y1, Image<unsigned char> image2, int x2, int y2, int halfWindowSize);

IndMatches vec_Matches;
features::PointFeatures feats0;
features::PointFeatures feats1;
Image<unsigned char> image0;
Image<unsigned char> image1;
int halfWindowSize = 3;

struct cmp{
    bool operator() ( IndMatch a, IndMatch b ){
        IndexT a_left = a._i, a_right = a._j;
        IndexT b_left = b._i, b_right = b._j;
        features::PointFeature a_left_feat = feats0[a_left];
        features::PointFeature a_right_feat = feats0[a_right];
        features::PointFeature b_left_feat = feats1[b_left];
        features::PointFeature b_right_feat = feats1[b_right];

        double a_zncc = zncc(image0, (int)a_left_feat.x(), (int)a_left_feat.y(), image1, (int)a_right_feat.x(), (int)a_right_feat.y(), halfWindowSize);
        double b_zncc = zncc(image0, (int)b_left_feat.x(), (int)b_left_feat.y(), image1, (int)b_right_feat.x(), (int)b_right_feat.y(), halfWindowSize);

        return a_zncc > b_zncc; }
};


int main() {
    const std::string sInputDir = stlplus::folder_up(string(THIS_SOURCE_DIR))
                                  + "imageData/SFMDir/";
    Hash_Map<IndexT, Image<unsigned char>> image_set;
    Hash_Map<IndexT, Image<RGBColor>> imageRGB_set;
    SfM_Data sfm_data;
    Mat3 K;
    Features_Provider feats_provider;
    Matches_Provider matches_provider;
    vector<std::string> filesname;

    DIR *dir;
    struct dirent *ent;
    int view_num = 0;


    if ((dir = opendir (sInputDir.c_str())) != NULL) {
        // read all the files within directory
        int i = 0;
        std::string image_path;
        Image<unsigned char> image;
        Image<RGBColor> imageRGB;
        while ((ent = readdir (dir)) != NULL) {
            if (stlplus::extension_part(ent->d_name) != "jpeg" && stlplus::extension_part(ent->d_name) != "jpg") {
                continue;
            }
            image_path = sInputDir + ent->d_name;
            filesname.push_back(image_path);
            ReadImage(image_path.c_str(), &image);
            ReadImage(image_path.c_str(), &imageRGB);
            image_set[i] = image;
            imageRGB_set[i] = imageRGB;
            i++;
            printf ("%s\n", image_path.c_str());
        }
        closedir (dir);
        view_num = i;
        if (readIntrinsic(image_path, K) < 0)
        {
            std::cerr << "Cannot read intrinsic parameters." << std::endl;
            return EXIT_FAILURE;
        }
    } else {
        // could not open directory
        perror ("");
        return EXIT_FAILURE;
    }


    using namespace openMVG::features;
    std::unique_ptr<Image_describer> image_describer(new SIFT_Image_describer);
    std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
    std::map<IndexT, std::map<int, int>> feature_to_track;


    sfm_data.views[0].reset(new View("", 0, 0, 0, image_set[0].Width(), image_set[0].Height()));
    sfm_data.intrinsics[0].reset(new Pinhole_Intrinsic(image_set[0].Width(), image_set[0].Height(), K(0, 0), K(0, 2), K(1, 2)));
    sfm_data.poses[sfm_data.views[0]->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    image_describer->Describe(image_set[0], regions_perImage[0]);

    /********** match adjacent two images **********/
    for(int k = 0; k < view_num - 1; k++) {
        const SIFT_Regions *regionsL = dynamic_cast<SIFT_Regions *>(regions_perImage.at(k).get());
        const PointFeatures featsL = regions_perImage.at(k)->GetRegionsPositions();
        if(k == 0) {
            feats_provider.feats_per_view[k] = featsL;
        }

        for (int i = k + 1; i < k + 2; i++) {
            cout << "************************************************" << endl
                 << "************* match " << k << " and " << i << "**************" << endl
                 << "************************************************"
                    << endl << endl;

            sfm_data.views[i].reset(new View("", i, 0, i, image_set[i].Width(), image_set[i].Height()));
            image_describer->Describe(image_set[i], regions_perImage[i]);

            const SIFT_Regions *regionsR = dynamic_cast<SIFT_Regions *>(regions_perImage.at(i).get());

            const PointFeatures featsR = regions_perImage.at(i)->GetRegionsPositions();
            feats_provider.feats_per_view[i] = featsR;

            /*
            // Draw features on two images
            {
                Image<unsigned char> concat;
                ConcatH(image_set[i - 1], image_set[i], concat);
                //-- Draw features :
                for (size_t m = 0; m < featsL.size(); ++m) {
                    const SIOPointFeature point = regionsL->Features()[m];
                    DrawCircle(point.x(), point.y(), point.scale(), 255, &concat);
                }
                for (size_t m = 0; m < featsR.size(); ++m) {
                    const SIOPointFeature point = regionsR->Features()[m];
                    DrawCircle(point.x() + image_set[i - 1].Width(), point.y(), point.scale(), 255, &concat);
                }
                std::string out_filename = "0_";
                out_filename.append(std::to_string(i)).append("_features.jpg");
                WriteImage(out_filename.c_str(), concat);
            }
             */

            std::vector<IndMatch> vec_PutativeMatches;
            //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
            {
                // Find corresponding points
                matching::DistanceRatioMatch(
                        0.8, matching::BRUTE_FORCE_L2,
                        *regions_perImage.at(k).get(),
                        *regions_perImage.at(i).get(),
                        vec_PutativeMatches);

                IndMatchDecorator<float> matchDeduplicator(
                        vec_PutativeMatches, featsL, featsR);
                matchDeduplicator.getDeduplicated(vec_PutativeMatches);

                matches_provider._pairWise_matches.insert(std::pair< Pair, IndMatches >(Pair(k, i), vec_PutativeMatches));

                std::cout
                << regions_perImage.at(k)->RegionCount() << " #Features on image " << k << std::endl
                << regions_perImage.at(i)->RegionCount() << " #Features on image " << i << std::endl
                << vec_PutativeMatches.size() << " #matches with Distance Ratio filter" << std::endl;

                // Draw correspondences after Nearest Neighbor ratio filter
                svgDrawer svgStream( image_set[k].Width() + image_set[i].Width(), max(image_set[k].Height(), image_set[i].Height()));
                svgStream.drawImage(filesname[k], image_set[k].Width(), image_set[k].Height());
                svgStream.drawImage(filesname[i], image_set[i].Width(), image_set[i].Height(), image_set[k].Width());
                for (size_t i = 0; i < vec_PutativeMatches.size(); ++i) {
                    //Get back linked feature, draw a circle and link them by a line
                    const SIOPointFeature L = regionsL->Features()[vec_PutativeMatches[i]._i];
                    const SIOPointFeature R = regionsR->Features()[vec_PutativeMatches[i]._j];
                    svgStream.drawLine(L.x(), L.y(), R.x()+image_set[k].Width(), R.y(), svgStyle().stroke("green", 2.0));
                    svgStream.drawCircle(L.x(), L.y(), L.scale(), svgStyle().stroke("yellow", 2.0));
                    svgStream.drawCircle(R.x()+image_set[k].Width(), R.y(), R.scale(), svgStyle().stroke("yellow", 2.0));
                }
                const std::string out_filename = "./Custom_SFM/" + std::to_string(k) + "_" + std::to_string(i) + "_siftMatches.svg";
                std::ofstream svgFile( out_filename.c_str() );
                svgFile << svgStream.closeSvgFile().str();
                svgFile.close();

            }

            //A. prepare the corresponding putatives points
            Mat xL(2, vec_PutativeMatches.size());
            Mat xR(2, vec_PutativeMatches.size());
            for (size_t m = 0; m < vec_PutativeMatches.size(); ++m) {
                const PointFeature &imaL = featsL[vec_PutativeMatches[m]._i];
                const PointFeature &imaR = featsR[vec_PutativeMatches[m]._j];
                xL.col(m) = imaL.coords().cast<double>();
                xR.col(m) = imaR.coords().cast<double>();
                //cout << m << ": " << endl << xL.col(m) << endl << xR.col(m) << endl;
            }

            //B. Compute the relative pose thanks to a essential matrix estimation
            std::pair<size_t, size_t> size_imaL(image_set[k].Width(), image_set[k].Height());
            std::pair<size_t, size_t> size_imaR(image_set[i].Width(), image_set[i].Height());
            RelativePose_Info relativePose_info;
            /*cout << "K: " << endl << K << endl;
            cout << "size:" << endl
            << image_set[k].Width() << " x " << image_set[k].Height() << endl
            << image_set[i].Width() << " x " << image_set[i].Height() << endl;*/
            if (!robustRelativePose(K, K, xL, xR, relativePose_info, size_imaL, size_imaR, 256)) {
                std::cerr << " /!\\ Robust relative pose estimation failure."
                << std::endl;
                return EXIT_FAILURE;
            }

            std::cout << "\nFound an Essential matrix:\n"
            << "\tprecision: " << relativePose_info.found_residual_precision << " pixels\n"
            << "\t#inliers: " << relativePose_info.vec_inliers.size() << "\n"
            << "\t#matches: " << vec_PutativeMatches.size()
            << std::endl;

            // Setup poses camera data
            sfm_data.poses[sfm_data.views[i]->id_pose] = relativePose_info.relativePose * sfm_data.poses[sfm_data.views[k]->id_pose];

            cout << "pose number " << sfm_data.views[k]->id_pose << ": \n"
            << "-- Rotation|Translation|Center matrices: --" << "\n"
            << sfm_data.poses[sfm_data.views[k]->id_pose].rotation() << "\n\n"
            << sfm_data.poses[sfm_data.views[k]->id_pose].translation() << "\n\n"
            << sfm_data.poses[sfm_data.views[k]->id_pose].center() << "\n"
            << std::endl << std::endl;

            cout << "relative pose from " << sfm_data.views[k]->id_pose << " to " << sfm_data.views[i]->id_pose << ": \n"
            << "-- Rotation|Translation|Center matrices: --" << "\n"
            << relativePose_info.relativePose.rotation() << "\n\n"
            << relativePose_info.relativePose.translation() << "\n\n"
            << relativePose_info.relativePose.center() << "\n"
            << std::endl << std::endl;

            cout << "pose number " << sfm_data.views[i]->id_pose << ": \n"
            << "-- Rotation|Translation|Center matrices: --" << "\n"
            << sfm_data.poses[sfm_data.views[i]->id_pose].rotation() << "\n\n"
            << sfm_data.poses[sfm_data.views[i]->id_pose].translation() << "\n\n"
            << sfm_data.poses[sfm_data.views[i]->id_pose].center() << "\n"
            << std::endl;


            //C. Triangulate and check valid points
            // invalid points that do not respect cheirality are discarded (removed
            //  from the list of inliers).
            const Pose3 pose0 = sfm_data.poses[sfm_data.views[k]->id_pose];
            const Pose3 pose1 = sfm_data.poses[sfm_data.views[i]->id_pose];

            // Init structure by inlier triangulation
            //get projection matrix
            const Mat34 P1 = sfm_data.intrinsics[sfm_data.views[k]->id_intrinsic]->get_projective_equivalent(pose0);
            const Mat34 P2 = sfm_data.intrinsics[sfm_data.views[i]->id_intrinsic]->get_projective_equivalent(pose1);
            Landmarks &landmarks = sfm_data.structure;
            for (size_t j = 0; j < relativePose_info.vec_inliers.size(); ++j) {
                const SIOPointFeature &LL = regionsL->Features()[vec_PutativeMatches[relativePose_info.vec_inliers[j]]._i];
                const SIOPointFeature &RR = regionsR->Features()[vec_PutativeMatches[relativePose_info.vec_inliers[j]]._j];
                std::map<IndexT, std::map<int, int>>::iterator ite;
                ite = feature_to_track.find(k);
                if (ite == feature_to_track.end()) {
                    feature_to_track[k] = std::map<int, int>();
                }
                ite = feature_to_track.find(i);
                if (ite == feature_to_track.end()) {
                    feature_to_track[i] = std::map<int, int>();
                }
                std::map<int, int>::iterator ite0;
                ite0 = feature_to_track[k].find(vec_PutativeMatches[relativePose_info.vec_inliers[j]]._i);
                if (ite0 != feature_to_track[k].end()) {
                    int track_id = feature_to_track[k][vec_PutativeMatches[relativePose_info.vec_inliers[j]]._i];
                    landmarks[track_id].obs[sfm_data.views[i]->id_view]
                            = Observation(RR.coords().cast<double>(),
                                          vec_PutativeMatches[relativePose_info.vec_inliers[j]]._j);
                    feature_to_track[i][vec_PutativeMatches[relativePose_info.vec_inliers[j]]._j] = track_id;
                }
                else {
                    // Point triangulation
                    int new_track_id = landmarks.size();
                    Vec3 X;
                    TriangulateDLT(P1, LL.coords().cast<double>(), P2, RR.coords().cast<double>(), &X);
                    // Reject point that is behind the camera
                    if (pose0.depth(X) < 0 && pose1.depth(X) < 0)
                        continue;
                    // Add a new landmark (3D point with it's 2d observations)
                    landmarks[new_track_id].obs[sfm_data.views[k]->id_view] = Observation(
                            LL.coords().cast<double>(), vec_PutativeMatches[relativePose_info.vec_inliers[j]]._i);
                    landmarks[new_track_id].obs[sfm_data.views[i]->id_view] = Observation(RR.coords().cast<double>(),
                                                                                          vec_PutativeMatches[relativePose_info.vec_inliers[j]]._j);
                    landmarks[new_track_id].X = X;
                    RGBColor RGBL = imageRGB_set[k]((int)LL.coords()[1], (int)LL.coords()[0]);
                    RGBColor RGBR = imageRGB_set[i]((int)LL.coords()[1], (int)LL.coords()[0]);

                    landmarks[new_track_id].RGB = Vec3((RGBL.r() + RGBR.r()) / 2, (RGBL.g() + RGBR.g()) / 2, (RGBL.b() + RGBR.b()) / 2);
                            //(imageRGB_set[k]((int)LL.coords()[1], (int)LL.coords()[0]) + imageRGB_set[i]((int)RR.coords()[1], (int)RR.coords()[0])) / 2;
                    feature_to_track[k][vec_PutativeMatches[relativePose_info.vec_inliers[j]]._i] = new_track_id;
                    feature_to_track[i][vec_PutativeMatches[relativePose_info.vec_inliers[j]]._j] = new_track_id;
                }
            }

            //D. Perform Bundle Adjustment of the scene

            Bundle_Adjustment_Ceres bundle_adjustment_obj;
            bundle_adjustment_obj.Adjust(sfm_data);

            std::string output = "./Custom_SFM/sparse_" + std::to_string(k) + "_" + std::to_string(i) + ".ply";
            Save(sfm_data, output, ESfM_Data(ALL));
        }
    }

    //perform dense match
    vec_Matches = matches_provider._pairWise_matches.at(Pair(0, 1));
    feats0 = regions_perImage.at(0)->GetRegionsPositions();
    feats1 = regions_perImage.at(1)->GetRegionsPositions();
    image0 = image_set.at(0);
    image1 = image_set.at(1);

    Mat matched0 = Mat::Zero(image0.Height(), image0.Width());
    Mat matched1 = Mat::Zero(image1.Height(), image1.Width());

    cout << "matched0" << ", [" << matched0.rows() << ", " << matched0.cols() << "]" << endl;
    cout << "15: [" << (int)feats0[vec_Matches[15]._i].x() << ", " << (int)feats0[vec_Matches[15]._i].y() << "], " <<
    "[" << (int)feats1[vec_Matches[15]._j].x() << ", " << (int)feats1[vec_Matches[15]._j].y() << "]\n";

    std::priority_queue<IndMatch, vector<IndMatch>, cmp> pq_matches;
    int cnt = 0;
    int good_matches = 0;
    for(IndMatch match : vec_Matches) {
        pq_matches.push(match);
        cout << cnt++ <<
                ", [" << (int)feats0[match._i].x() << ", " << (int)feats0[match._i].y() << "], " <<
                "[" << (int)feats1[match._j].x() << ", " << (int)feats1[match._j].y() << "], ";
        matched0((int)feats0[match._i].y(), (int)feats0[match._i].x()) = 1;
        matched1((int)feats1[match._j].y(), (int)feats1[match._j].x()) = 1;
        double cur_zncc = zncc(image0, (int)feats0[match._i].x(), (int)feats0[match._i].y(), image1, (int)feats1[match._j].x(), (int)feats1[match._j].y(), halfWindowSize);
        if(cur_zncc > 0.5) {
            good_matches ++;
        }
        cout <<"correlation between matched points: " << cur_zncc << endl;
    }
    cout << good_matches <<" good matches in " << vec_Matches.size() << " total matches." << endl;

    cout << "zncc test: " << zncc(image0, 100, 100, image0, 100, 100, halfWindowSize) << endl;

/*    while(!pq_matches.empty()) {
        IndMatch cur_match = pq_matches.top();
        pq_matches.pop();

        int x0 = (int)feats0[cur_match._i].x(), y0 = (int)feats0[cur_match._i].y(),
                x1 = (int)feats1[cur_match._j].x(), y1 = (int)feats1[cur_match._j].y();

        int xmin0 = std::max(halfWindowSize + 1, x0 - halfWindowSize),
                xmax0 = std::min(image0.Width() - halfWindowSize, x0 + halfWindowSize + 1),
                ymin0 = std::max(halfWindowSize + 1, y0 - halfWindowSize),
                ymax0 = std::min(image0.Height() - halfWindowSize, y0 + halfWindowSize + 1);

        int xmin1 = std::max(halfWindowSize + 1, x1 - halfWindowSize),
                xmax1 = std::min(image1.Width() - halfWindowSize, x1 + halfWindowSize + 1),
                ymin1 = std::max(halfWindowSize + 1, y1 - halfWindowSize),
                ymax1 = std::min(image1.Height() - halfWindowSize, y1 + halfWindowSize + 1);

        for(int xx0 = xmin0; xx0 <= xmax0; xx0++) {
            for(int yy0 = ymin0; yy0 <= ymax0; yy0++) {
                if(matched0(xx0, yy0)) {continue;}
                int xx = xx0 - x0 + x1;
                int yy = yy0 - y0 + y1;
                for(int yy1 = std::max(ymin1, yy - 1); yy1 <= std::min(ymax1, yy + 2); yy1++) {
                    for(int xx1 = std::max(xmin1, xx - 1); yy1 <= std::min(xmax1, xx + 2); xx1++) {
                        if(!matched1(xx1, yy1)) {
                            double correlation = zncc(image0, xx0, yy0, image1, xx1, yy1, halfWindowSize);
                            if(correlation > 0.8) {
                                PointFeature point0(xx0, yy0), point1(xx1, yy1);
                                feats0.push_back(point0);
                                feats1.push_back(point1);
                                matched0(xx0, yy0) = 1;
                                matched1(xx1, yy1) = 1;

                                IndMatch new_match(feats0.size() - 1, feats1.size() - 1);
                                pq_matches.push(new_match);

                            }
                        }
                    }
                }

            }
        }

    } */





/*
    cout << "************************************************" << endl
    << "**************** Use Sequential SFM Pipeline **************" << endl
    << "************************************************"
    << endl << endl;

    SfM_Data sfm_data_2 = sfm_data;
    sfm_data_2.poses.clear();
    sfm_data_2.structure.clear();

    SequentialSfMReconstructionEngine sfmEngine(
            sfm_data_2,
            "./Sequential_SFM/",
            stlplus::create_filespec("./Sequential_SFM/", "Sequential_Reconstruction_Report.html"));


    // Configure data provider (Features and Matches)
    sfmEngine.SetFeaturesProvider(&feats_provider);
    sfmEngine.SetMatchesProvider(&matches_provider);

    // Set an initial pair
    sfmEngine.setInitialPair(Pair(0,1));

    // Configure reconstruction parameters
    sfmEngine.Set_bFixedIntrinsics(true);

    sfmEngine.Process();


    cout << "************************************************" << endl
    << "**************** Use Global SFM Pipeline **************" << endl
    << "************************************************"
    << endl << endl;

    SfM_Data sfm_data_3 = sfm_data;
    sfm_data_3.poses.clear();
    sfm_data_3.structure.clear();

    GlobalSfMReconstructionEngine_RelativeMotions global_sfmEngine(
            sfm_data_3,
            "./Global_SFM/",
            stlplus::create_filespec("./Global_SFM/", "Global_Reconstruction_Report.html"));


    // Configure data provider (Features and Matches)
    global_sfmEngine.SetFeaturesProvider(&feats_provider);
    global_sfmEngine.SetMatchesProvider(&matches_provider);

    // Configure reconstruction parameters
    global_sfmEngine.Set_bFixedIntrinsics(true);

    // Configure motion averaging method
    global_sfmEngine.SetRotationAveragingMethod(ROTATION_AVERAGING_L2);
    global_sfmEngine.SetTranslationAveragingMethod(TRANSLATION_AVERAGING_L1);

    global_sfmEngine.Process();
    */

}


int readIntrinsic(const std::string & fileName, Mat3 & K)
{
    // Load the K matrix
    FILE *fp = fopen(fileName.c_str(), "rb");
    if(!fp) {
        cout << "Cannot open " << fileName << endl;
        return -1;
    }
    fseek(fp, 0, SEEK_END);
    unsigned long fsize = ftell(fp);
    rewind(fp);
    unsigned char *buf = new unsigned char[fsize];
    if (fread(buf, 1, fsize, fp) != fsize) {
        printf("Can't read file.\n");
        delete[] buf;
        return -2;
    }
    fclose(fp);

    EXIFInfo result;
    int code = result.parseFrom(buf, fsize);
    delete[] buf;
    if (code) {
        printf("Error parsing EXIF: code %d\n", code);
        return -3;
    }


    double CCD_width = 0;
    double f = 0;
    //take care of apple
    if (strcmp(result.Make.c_str(), "Apple") == 0 && result.FocalLengthIn35mm == 35) {
        CCD_width = 36;
        int width = (result.ImageWidth > result.ImageHeight) ? result.ImageWidth : result.ImageHeight;
        f = result.FocalLengthIn35mm * width / CCD_width;
    }
    K(0, 0) = f;
    K(1, 1) = f;
    K(2, 2) = 1;
    K(0, 2) = result.ImageWidth / 2;
    K(1, 2) = result.ImageHeight / 2;

    cout << "intrinsic matrix:" << endl << K << endl;
    return 0;

}

double zncc(Image<unsigned char> image1, int x1, int y1, Image<unsigned char> image2, int x2, int y2, int halfWindowSize) {
/*    int xmin1 = std::max(0, x1 - halfWindowSize);
    int xmax1 = std::min(image1.Width() - 1, x1 + halfWindowSize);
    int xmin2 = std::max(0, x2 - halfWindowSize);
    int xmax2 = std::min(image2.Width() - 1, x2 + halfWindowSize);
    int ymin1 = std::max(0, y1 - halfWindowSize);
    int ymax1 = std::min(image1.Height() - 1, y1 + halfWindowSize);
    int ymin2 = std::max(0, y2 - halfWindowSize);
    int ymax2 = std::min(image2.Height() - 1, y2 + halfWindowSize);
    int xmin = std::max(xmin1, xmin2);*/

    if(x1 - halfWindowSize < 0 || x1 + halfWindowSize > image1.Width() - 1) {
        return -2;
    }
    if(x2 - halfWindowSize < 0 || x2 + halfWindowSize > image2.Width() - 1) {
        return -2;
    }
    if(y1 - halfWindowSize < 0 || y1 + halfWindowSize > image1.Height() - 1) {
        return -2;
    }
    if(y2 - halfWindowSize < 0 || y2 + halfWindowSize > image2.Height() - 1) {
        return -2;
    }

    double SI = 0, SJ = 0, SII = 0, SJJ = 0, SIJ = 0;
    for(int i = -halfWindowSize; i <= halfWindowSize; i++) {
        for(int j = -halfWindowSize; j <= halfWindowSize; j++) {
            int xx1 = x1 + i, yy1 = y1 + j;
            int xx2 = x2 + i, yy2 = y2 + j;
            if(xx1 < 0 || xx1 > image1.Width() || yy1 < 0 || yy1 > image1.Height()) {
                return -3;
            }
            if(xx2 < 0 || xx2 > image2.Width() || yy2 < 0 || yy2 > image2.Height()) {
                return -3;
            }
            SI += image1(yy1, xx1) / 256.0;
            SJ += image2(yy2, xx2) / 256.0;
            SII += std::pow(image1(yy1, xx1) / 256.0, 2);
            SJJ += std::pow(image2(yy2, xx2) / 256.0, 2);
            SIJ += (image1(yy1, xx1) / 256.0) * (image2(yy2, xx2) / 256.0);
        }
    }

    int N = std::pow(2 * halfWindowSize + 1, 2);
    double ret = (N * SIJ - SI * SJ) / std::sqrt((N * SII - SI * SI) * (N * SJJ - SJ * SJ));

    return ret;
}



