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

    vec3 relativeViewPosition = vec4(viewTransform * vec4(0.0F, 1.0F)).xyz;
    vec3 relativePosition = vec4(viewTransform * vec4(vertCoordinates_in, height, 1.0F)).xyz;
    vec3 relativeLightPosition = vec4(viewTransform * lightPosition).xyz;

    vec3 N = normalize(normalTransform * vertNormals_in);
    vec3 L = normalize(relativeLightPosition - relativePosition);
    vec3 R = 2*dot(L, N)*N - L;
    vec3 V = normalize(relativeViewPosition - relativePosition);

    shading = material[0] + material[1]*max(dot(N, L), 0.0F) + material[2]*pow(max(dot(R, V), 0.0F), material[3]);
    heightChange = length(vertNormals_in.xy);
}
