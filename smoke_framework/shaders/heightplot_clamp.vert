#version 330 core
// Height plot clamp vertex shader

layout (location = 0) in vec2 vertCoordinates_in;
layout (location = 1) in float height;
layout (location = 2) in float value_in;
layout (location = 3) in vec3 vertNormals_in;

// height value
out float value;
// depicts the brightness of the color
out float shading;
// magnitude of the gradient
out float heightChange;

uniform float clampMin;
uniform float clampMax;
uniform float transferK;

uniform mat4 projectionTransform;
uniform mat4 viewTransform;
uniform mat3 normalTransform;

uniform vec4 material; // Contains 4 floats, in order: k_a, k_d, k_s, alpha.
uniform vec3 lightPosition;

void main()
{
    // Clamp height
    value = clamp(height, clampMin, clampMax);

    gl_Position = viewTransform * projectionTransform * vec4(vertCoordinates_in, height, 1.0F);

    vec3 viewPosition = vec3(0.0F);
    vec3 relativeVertexPosition = vec3(viewTransform * vec4(vertCoordinates_in, height, 1.0F));
    vec3 relativeLightPosition = vec3(viewTransform * vec4(lightPosition, 1.0F));

    vec3 N = normalize(normalTransform * vertNormals_in);
    vec3 L = normalize(relativeLightPosition - relativeVertexPosition);
    vec3 R = 2*dot(L, N)*N - L;
    vec3 V = normalize(viewPosition - relativeVertexPosition);

    shading  = material[0];
    shading += material[1] * max(dot(N, L), 0.0F);
    shading += material[2] * pow(max(dot(R, V), 0.0F), material[3]);

    heightChange = length(vertNormals_in.xy);
}
