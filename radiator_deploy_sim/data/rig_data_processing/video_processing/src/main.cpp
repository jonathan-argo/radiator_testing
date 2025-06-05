#include <opencv2/opencv.hpp>  
#include <opencv2/aruco.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int main() {

    cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_250);

    #ifdef GENERATE_MARKERS
        cv::Mat markerImage;

        for (int i = 0; i < 4; ++i) {
            int markerId = 0 + i;
            cv::Mat markerImage;
            cv::aruco::generateImageMarker(dictionary, markerId, 200, markerImage, 1);

            std::string filename = "../../fiducials/marker_" + std::to_string(markerId) + ".png";
            cv::imwrite(filename, markerImage);
        }
    #endif

    #ifdef SQUARE_VIDEO
        cv::VideoCapture cap("/home/jonathanvollrath/radiator_testing/radiator_deploy_sim/data/test_video.mp4");
        if (!cap.isOpened()) {
            std::cerr << "Could not open video file.\n";
            return -1;
        }

        cv::VideoWriter writer;
        int codec = cv::VideoWriter::fourcc('m', 'p', '4', 'v');
        double fps = cap.get(cv::CAP_PROP_FPS);
        cv::Size outputSize(500, 500); // this will need to be updated for the size of the video file that the visualization spits out

        writer.open("/home/jonathanvollrath/radiator_testing/radiator_deploy_sim/data/test_output.mp4", codec, fps, outputSize, true);
        if (!writer.isOpened()) {
            std::cerr << "Failed to open video writer.\n";
            return -1;
        }

        // ArUco dictionary and detector
        cv::aruco::DetectorParameters detectorParams;
        cv::aruco::ArucoDetector detector(dictionary, detectorParams);

        // Target size of the rectified image
        int outputWidth = 500, outputHeight = 500;

        while (true) {
            cv::Mat frame;
            cap >> frame;
            if (frame.empty()) break;

            std::vector<int> ids;
            std::vector<std::vector<cv::Point2f>> corners, rejected;
            detector.detectMarkers(frame, corners, ids, rejected);

            if (ids.size() >= 4) {
                std::map<int, cv::Point2f> markerCenters;

                for (size_t i = 0; i < ids.size(); ++i) {
                    if (ids[i] >= 0 && ids[i] <= 3) { // Only use markers 0,1,2,3
                        cv::Point2f center(0, 0);
                        for (const auto& pt : corners[i])
                            center += pt;
                        center *= 0.25f;
                        markerCenters[ids[i]] = center;
                    }
                }

                if (markerCenters.size() == 4) {
                    std::vector<cv::Point2f> src = {
                        markerCenters[0],
                        markerCenters[1],
                        markerCenters[3],
                        markerCenters[2]
                    };
                    std::vector<cv::Point2f> dst = {
                        {0.0f, 0.0f},
                        {(float)outputWidth, 0.0f},
                        {(float)outputWidth, (float)outputHeight},
                        {0.0f, (float)outputHeight}
                    };

                    cv::Mat H = cv::getPerspectiveTransform(src, dst);
                    cv::Mat warped;
                    cv::warpPerspective(frame, warped, H, cv::Size(outputWidth, outputHeight));

                    writer.write(warped); 
                    cv::imshow("Warped", warped);
                }
            }

            cv::imshow("Original", frame);
            char key = (char) cv::waitKey(30);
            if (key == 27) break; // ESC to quit
        }
    #endif

    return 0;
}
