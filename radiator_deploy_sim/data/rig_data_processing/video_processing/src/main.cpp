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

    cv::VideoCapture inputVideo;
    inputVideo.open("/home/jonathanvollrath/radiator_testing/radiator_deploy_sim/data/test_video.mp4");

    if (!inputVideo.isOpened()) {
        std::cerr << "Could not open the video file." << std::endl;
        return -1;
    }

    cv::aruco::DetectorParameters detectorParams = cv::aruco::DetectorParameters();
    cv::aruco::ArucoDetector detector(dictionary, detectorParams);

    const int waitTime = 30; // milliseconds between frames

    while (inputVideo.grab()) {
        cv::Mat image, imageCopy;
        inputVideo.retrieve(image);
        image.copyTo(imageCopy);

        std::vector<int> ids;
        std::vector<std::vector<cv::Point2f>> corners, rejected;
        detector.detectMarkers(image, corners, ids, rejected);

        if (!ids.empty())
            cv::aruco::drawDetectedMarkers(imageCopy, corners, ids);

        cv::imshow("out", imageCopy);
        char key = (char) cv::waitKey(waitTime);
        if (key == 27) // ESC
            break;
    }
    return 0;
}
