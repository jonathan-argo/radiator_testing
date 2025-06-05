#include <opencv2/aruco.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>

int main() {

    cv::Mat markerImage;
    cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_250);

    for (int i = 0; i < 4; ++i) {
        int markerId = 23 + i;
        cv::Mat markerImage;
        cv::aruco::generateImageMarker(dictionary, markerId, 200, markerImage, 1);

        std::string filename = "../../fiducials/marker_" + std::to_string(markerId) + ".png";
        cv::imwrite(filename, markerImage);
    }

    return 0;
}
