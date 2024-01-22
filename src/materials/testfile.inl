#include <cmath>

// Function to calculate the Fresnel term F_m using the Schlick approximation
float calculateFresnel(float baseColor, const glm::vec3& h, const glm::vec3& omegaOut) {
    return baseColor + (1.0f - baseColor) * pow(1.0f - glm::dot(h, omegaOut), 5.0f);
}

// Function to calculate the normal distribution function D_m (GGX)
float calculateNormalDistribution(float alphaX, float alphaY, const glm::vec3& h) {
    float hxl = h.x / alphaX;
    float hyl = h.y / alphaY;
    float hzl = h.z;

    return 1.0f / (M_PI * alphaX * alphaY * pow(hxl * hxl + hyl * hyl + hzl * hzl, 2.0f));
}

// Function to calculate the aspect, alpha_x, and alpha_y
void calculateAspectAndAlphas(float anisotropic, float roughness, float& aspect, float& alphaX, float& alphaY) {
    aspect = sqrt(1.0f - 0.9f * anisotropic);
    float alphaMin = 0.0001f;

    alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
}

// Function to calculate the average occlusion factor G_m using the Smith model
float calculateOcclusionFactor(const glm::vec3& omegaIn, const glm::vec3& omegaOut, float alphaX, float alphaY) {
    auto lambda = [&](const glm::vec3& omegaL) {
        return sqrt(1.0f + pow((omegaL.x * alphaX), 2.0f) + pow((omegaL.y * alphaY), 2.0f) / pow(omegaL.z, 2.0f)) - 1.0f;
    };

    auto G = [&](const glm::vec3& omega) {
        return 1.0f / (1.0f + lambda(omega) / 2.0f);
    };

    return G(omegaIn) * G(omegaOut);
}

// Function to calculate the specular reflection for metal
float calculateMetalSpecularReflection(float baseColor, float alphaX, float alphaY, const glm::vec3& n, const glm::vec3& omegaIn) {
    glm::vec3 h = glm::normalize(omegaIn + n);
    float fresnel = calculateFresnel(baseColor, h, omegaIn);
    float normalDistribution = calculateNormalDistribution(alphaX, alphaY, h);
    float occlusionFactor = calculateOcclusionFactor(omegaIn, omegaIn, alphaX, alphaY);

    return (fresnel * normalDistribution * occlusionFactor) / (4.0f * std::abs(glm::dot(n, omegaIn)));
}
