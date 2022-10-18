#include "visualization.h"

#include "constants.h"
#include "mainwindow.h"

#include <QDebug>
#include <QVector3D>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

Visualization::Visualization(QWidget *parent) : QOpenGLWidget(parent)
{
    qDebug() << "Visualization constructor";

    // Start the simulation loop.
    m_timer.start(17); // Each frame takes 17ms, making the simulation run at approximately 60 FPS
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(do_one_simulation_step()));

    m_elapsedTimer.start();

    std::fill(m_scalarCube.begin(), m_scalarCube.end(), std::vector<float>(m_DIM * m_DIM, 0.0F));
}

Visualization::~Visualization()
{
    makeCurrent();

    qDebug() << "Visualization destructor";

    opengl_deleteObjects();
}

void Visualization::do_one_simulation_step()
{
    if (m_isRunning)
        m_simulation.do_one_simulation_step();

    update();
}

void Visualization::initializeGL() {
    qDebug() << ":: Initializing OpenGL";
    initializeOpenGLFunctions();

    connect(&m_debugLogger, SIGNAL(messageLogged(QOpenGLDebugMessage)),
            this, SLOT(onMessageLogged(QOpenGLDebugMessage)), Qt::DirectConnection);

    if (m_debugLogger.initialize())
    {
        qDebug() << ":: Logging initialized";
        m_debugLogger.startLogging( QOpenGLDebugLogger::SynchronousLogging );
        m_debugLogger.enableMessages();
    }

    {
        QString const glVersion{reinterpret_cast<const char*>(glGetString(GL_VERSION))};
        qDebug() << ":: Using OpenGL" << qPrintable(glVersion);
    }

    glClearColor(0.2F, 0.1F, 0.2F, 1.0F);

    // Retrieve default textures.
    auto const mainWindowPtr = qobject_cast<MainWindow*>(parent()->parent());
    std::vector<Color> const defaultScalarDataColorMap = mainWindowPtr->m_defaultScalarDataColorMap;
    std::vector<Color> const defaultVectorDataColorMap = mainWindowPtr->m_defaultVectorDataColorMap;

    opengl_generateObjects();
    opengl_createShaderPrograms();

    opengl_setupAllBuffers();

    opengl_loadScalarDataTexture(defaultScalarDataColorMap);
    opengl_loadVectorDataTexture(defaultVectorDataColorMap);
    opengl_loadLicTexture(std::vector<uint8_t>()); // Initially provide an empty texture

    opengl_rotateView();
}

void Visualization::paintGL()
{
    glBindVertexArray(0U);

    // Clear the screen before rendering
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // The height plot, LIC and volume rendering must be drawn by themselves.
    // The scalar data, isolines and vector data drawing can be combined.
    if (m_drawHeightplot)
    {
        opengl_drawHeightplot();
        return;
    }

    if (m_drawLIC)
    {
        opengl_drawLic();
        return;
    }

    if (m_drawVolumeRendering)
    {
        opengl_drawVolumeRendering();
        return;
    }

    if (m_drawScalarData)
        drawScalarData();

    if (m_drawIsolines)
    {
        m_shaderProgramIsolines.bind();
        glUniformMatrix4fv(m_uniformLocationIsolines_projection, 1, GL_FALSE, m_projectionTransformationMatrix.data());
        glUniform3fv(m_uniformLocationIsolines_color, 1, &m_isolineColor[0]);
        opengl_drawIsolines();
    }

    if (m_drawVectorData)
    {
        m_shaderProgramVectorData.bind();
        glUniformMatrix4fv(m_uniformLocationProjectionColorMapInstanced, 1, GL_FALSE, m_projectionTransformationMatrix.data());
        glUniform1i(m_uniformLocationTextureColorMapInstanced, 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_1D, m_vectorDataTextureLocation);
        drawGlyphs();
    }
}

void Visualization::resizeGL(int const width, int const height)
{
    m_cellWidth  = static_cast<float>(width) / static_cast<float>(m_DIM + 1U);
    m_cellHeight = static_cast<float>(height) / static_cast<float>(m_DIM + 1U);

    m_projectionTransformationMatrix.setToIdentity();

    switch (m_projectionType)
    {
        case ProjectionType::Orthographic:
            m_projectionTransformationMatrix.ortho(0.0F, width,
                                                   0.0F, height,
                                                   -50.0F, 50.0F);
        break;

        case ProjectionType::Perspective:
            m_projectionTransformationMatrix.perspective(30.0F, 1.0F, 0.2F, 1300.0F);
            m_projectionTransformationMatrix.lookAt({300.0F, 300.0F, 1200.0F}, {300.0F, 300.0F, 0.0F}, {0.0F, 1.0F, 0.0F});
        break;
    }

    m_glyphCellWidth = static_cast<float>(width) / static_cast<float>(m_numberOfGlyphsX + 1U);
    m_glyphCellHeight = static_cast<float>(height) / static_cast<float>(m_numberOfGlyphsY + 1U);

    opengl_updateScalarPoints();
    opengl_updateLicPoints();
}

void Visualization::drawGlyphs()
{
    std::vector<float> vectorMagnitude;
    std::vector<float> vectorDirectionX;
    std::vector<float> vectorDirectionY;
    switch (m_currentVectorDataType)
    {
        case VectorDataType::Velocity:
            vectorMagnitude = m_simulation.velocityMagnitudeInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
            vectorDirectionX = m_simulation.velocityXInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
            vectorDirectionY = m_simulation.velocityYInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
        break;

        case VectorDataType::ForceField:
            vectorMagnitude = m_simulation.forceFieldMagnitudeInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
            vectorDirectionX = m_simulation.forceFieldXInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
            vectorDirectionY = m_simulation.forceFieldYInterpolated(m_numberOfGlyphsX, m_numberOfGlyphsY);
        break;
    }

    // Scale the magnitudes to where these become visible.
    std::transform(vectorMagnitude.begin(), vectorMagnitude.end(), vectorMagnitude.begin(),
                   std::bind(std::multiplies<>(), std::placeholders::_1, m_vectorDataMagnifier));

    if (m_sendMinMaxToUI && !vectorMagnitude.empty())
    {
        auto const currentMinMaxIt = std::minmax_element(vectorMagnitude.cbegin(), vectorMagnitude.cend());
        QVector2D const currentMinMax{*currentMinMaxIt.first, *currentMinMaxIt.second};

        // Send values to GUI.
        auto const mainWindowPtr = qobject_cast<MainWindow*>(parent()->parent());
        Q_ASSERT(mainWindowPtr != nullptr);
        mainWindowPtr->setVectorDataMin(currentMinMax.x());
        mainWindowPtr->setVectorDataMax(currentMinMax.y());
    }

    size_t const numberOfInstances = m_numberOfGlyphsX * m_numberOfGlyphsY;

    // Create model transformation matrices
    std::vector<float> modelTransformationMatrices;
    /* Fill the container modelTransformationMatrices here...
     * Use the following variables:
     * modelTransformationMatrix: This vector should contain the result.
     * m_DIM: The grid dimensions of the simulation.
     * m_cellWidth, m_cellHeight: A cell, made up of 4 simulation grid points, has the size m_cellWidth * m_cellHeight for the visualization.
     * m_numberOfGlyphsX (horizontal)
     * m_numberOfGlyphsY (vertical)
     * m_glyphCellWidth, m_glyphCellHeight: Use these as a small offset from the border of the OpenGL window.
     *                                      Having this little border around the visualization prevents a left-pointing arrow
     *                                      on the left side of the window to become invisible and thus convey no information.
     * vectorDirectionX, vectorDirectionY: To which direction the glyph should point. Row-major, size given by the m_numberOfGlyphs*.
     * vectorMagnitude: Use this value to scale the glyphs. I.e. higher values are visualized using larger glyphs. Row-major, size given by the m_numberOfGlyphs*.
     */
    modelTransformationMatrices = std::vector<float>(numberOfInstances * 16U, 0.0F); // Remove this placeholder initialization

    // TODO: This shouldn't be here, but otherwise re-binding an already bound Glyphs VAO causes glitches.
    glBindVertexArray(0);

    // Buffering section starts here.
    glBindVertexArray(m_vaoGlyphs);

    glBindBuffer(GL_ARRAY_BUFFER, m_vboValuesGlyphs);
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    static_cast<GLsizeiptr>(vectorMagnitude.size() * sizeof(float)),
                    vectorMagnitude.data());

    // Buffer model transformation matrices.
    glBindBuffer(GL_ARRAY_BUFFER, m_vboModelTransformationMatricesGlyphs);
    void * const dataPtr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    memcpy(dataPtr, modelTransformationMatrices.data(), modelTransformationMatrices.size() * sizeof(float));
    glUnmapBuffer(GL_ARRAY_BUFFER);

    if (m_currentGlyphType == Glyph::GlyphType::Hedgehog)
        glDrawElementsInstanced(GL_LINES,
                                static_cast<GLsizei>(m_glyphIndicesSize),
                                GL_UNSIGNED_SHORT,
                                reinterpret_cast<GLvoid*>(0),
                                static_cast<GLsizei>(numberOfInstances));
    else
        glDrawElementsInstanced(GL_TRIANGLE_STRIP,
                                static_cast<GLsizei>(m_glyphIndicesSize),
                                GL_UNSIGNED_SHORT,
                                reinterpret_cast<GLvoid*>(0),
                                static_cast<GLsizei>(numberOfInstances));
}


void Visualization::applyQuantization(std::vector<float> &scalarValues) const
{
    // Convert the floating point values to (8 bit) unsigned integers,
    // so that the data can be treated as an image.
    // The image's pixel values are in the range [0, 255].
    float const maxValue = *std::max_element(scalarValues.cbegin(), scalarValues.cend());
    std::vector<unsigned int> image;
    image.reserve(scalarValues.size());
    for (auto const x : scalarValues)
        image.push_back(static_cast<unsigned int>(std::lroundf(x / maxValue * 255.0F)));


    // Apply quantization to std::vector<unsigned int> image here.
    // The variable m_quantizationBits ('n' in the lecture slides) is set in the GUI and can be used here.
    // L needs to be set to the appropriate value and will be used to set the clamping range in the GUI.
    int block_size = 1;
    switch (m_quantizationBits) {
        case 1:
            block_size = 128;
            break;
        case 2:
            block_size = 64;
            break;
        case 4:
            block_size = 16;
            break;
        case 8:
            block_size = 1;
            break;
    }

    for (size_t i=0U; i<image.size(); i++) {
        image[i] = round(image[i] / block_size);
    }

    unsigned int const L = unsigned(256 / block_size) - 1;

    qDebug() << "Quantization not implemented";

    // Convert the image's data back to floating point values, so that it can be processed as usual.
    scalarValues = std::vector<float>{image.cbegin(), image.cend()};

    // Force the clamping range in the GUI to be [0, L].
    auto const mainWindowPtr = qobject_cast<MainWindow*>(parent()->parent());
    Q_ASSERT(mainWindowPtr != nullptr);
    mainWindowPtr->on_scalarDataMappingClampingMaxSlider_valueChanged(0);
    mainWindowPtr->on_scalarDataMappingClampingMaxSlider_valueChanged(100 * static_cast<int>(L));
}

std::vector<float> Visualization::applyKernel(std::vector<std::vector<float>> kernel, std::vector<float> image) const
{
    std::vector<float> newImage(m_DIM*m_DIM, 0.0F);
    int count;

    for (size_t x=0U; x<m_DIM; x++)
    {
        for (size_t y=0U; y<m_DIM; y++)
        {
            count = 0;

            if (x>0 && y>0) {
                newImage[m_DIM*x+y] += image[m_DIM*(x-1)+y-1] * kernel[0][0];
                count++;
            }
            if (x>0) {
                newImage[m_DIM*x+y] += image[m_DIM*(x-1)+y]   * kernel[0][1];
                count++;
            }
            if (x>0 && y<m_DIM-1) {
                newImage[m_DIM*x+y] += image[m_DIM*(x-1)+y+1] * kernel[0][2];
                count++;
            }
            if (y>0) {
                newImage[m_DIM*x+y] += image[m_DIM*x+y-1]     * kernel[1][0];
                count++;
            }
            if (1) {
                newImage[m_DIM*x+y] += image[m_DIM*x+y]       * kernel[1][1];
                count++;
            }
            if (y<m_DIM-1) {
                newImage[m_DIM*x+y] += image[m_DIM*x+y+1]     * kernel[1][2];
                count++;
            }
            if (x<m_DIM-1 && y>0) {
                newImage[m_DIM*x+y] += image[m_DIM*(x+1)+y-1] * kernel[2][0];
                count++;
            }
            if (x<m_DIM-1) {
                newImage[m_DIM*x+y] += image[m_DIM*(x+1)+y]   * kernel[2][1];
                count++;
            }
            if (x<m_DIM-1 && y<m_DIM-1) {
                newImage[m_DIM*x+y] += image[m_DIM*(x+1)+y+1] * kernel[2][2];
                count++;
            }

            newImage[m_DIM*x+y] /= count;
        }
    }

    return newImage;
}

void Visualization::applyGaussianBlur(std::vector<float> &scalarValues) const
{
    std::vector<float> NewImageX(m_DIM*m_DIM);
    std::vector<std::vector<float>> kernel = {{ 1, 2, 1 },{ 2, 4, 2 },{ 1, 2, 1 }};

    scalarValues = applyKernel(kernel, scalarValues);
}

void Visualization::applyGradients(std::vector<float> &scalarValues) const
{
    // Implement Gradient extraction here, applied on the values of the scalarValues container.
    // First, define a 3x3 Sobel kernels (for x and y directions).
    // (Use a C-style 2D array, a std::array of std::array's, or a std::vector of std::vectors)
    // Convolve the values of the scalarValues container with the Sobel kernels
    // Calculate the Gradient magnitude
    // Calculate the Gradient direction
    // Visualize the Gradient magnitude

    std::vector<std::vector<float>> xKernel = {{ 1, 0, -1 },{ 2, 0, -2 },{ 1, 0, -1 }};
    std::vector<std::vector<float>> yKernel = {{ 1, 2, 1 },{ 0, 0, 0 },{ -1, -2, -1 }};

    std::vector<float> xDir = applyKernel(xKernel, scalarValues);
    std::vector<float> yDir = applyKernel(yKernel, scalarValues);

    for (size_t i=1U; i<m_DIM*m_DIM; i++) {
        scalarValues[i] = sqrt(xDir[i]*xDir[i] + yDir[i]*yDir[i]);
    }
}

/* This function receives a *reference* to a std::vector<float>,
 * which acts as a pointer. Modifying scalarValues here will result
 * in the scalarValues passed to this function to be modified.
 * You may also define a new "std::vector<float> v" here, fill it with
 * new values and assign it to scalarValues to replace it,
 * e.g. "scalarValues = v;".
 *
 * m_sliceIdx contains the value set in the GUI.
 * m_DIM contains the current dimensions of the square (m_DIM * m_DIM).
 * m_slicingWindowSize contains the size of the window (here, also m_DIM).
 * m_slicingDirection contains the slicing direction set in the GUI and
 *    is already handled in a switch statement.
 */
void Visualization::applySlicing(std::vector<float> &scalarValues)
{
    // Shift space slices back in time, starting at end
    for (size_t t=m_DIM-1; t>0; t--) {
        m_scalarCube[t] = m_scalarCube[t-1];
    }

    // Add new space slice at begin
    m_scalarCube[0] = scalarValues;

    // Slice to render
    std::vector<float> slice;

    switch (m_slicingDirection)
    {
    case SlicingDirection::x:
        // xIdx is constant
        for (int t=0U; t<m_DIM; t++) {
            for (int y=0U; y<m_DIM; y++) {
                slice.push_back(m_scalarCube[t][m_sliceIdx*m_DIM + y]);
            }
        }
        break;

    case SlicingDirection::y:
        // yIdx is constant
        for (size_t x=0U; x<m_DIM; x++) {
            for (size_t t=0U; t<m_DIM; t++) {
                slice.push_back(m_scalarCube[t][x*m_DIM + m_sliceIdx]);
            }
        }
        break;

    case SlicingDirection::t:
        // t is constant
        for (size_t x=0U; x<m_DIM; x++) {
            for (size_t y=0U; y<m_DIM; y++) {
                slice.push_back(m_scalarCube[m_sliceIdx][x*m_DIM + y]);
            }
        }
        break;
    }

    scalarValues = slice;
}

void Visualization::applyPreprocessing(std::vector<float> &scalarValues)
{
    if (m_useQuantization)
        applyQuantization(scalarValues);

    if (m_useGaussianBlur)
        applyGaussianBlur(scalarValues);

    if (m_useGradients)
        applyGradients(scalarValues);

    if (m_useSlicing)
        applySlicing(scalarValues);
}

void Visualization::drawScalarData()
{
    std::vector<float> scalarValues;

    switch (m_currentScalarDataType)
    {
        case ScalarDataType::Density:
            scalarValues = m_simulation.density();
        break;

        case ScalarDataType::ForceFieldMagnitude:
            scalarValues = m_simulation.forceFieldMagnitude();
        break;

        case ScalarDataType::VelocityMagnitude:
            scalarValues = m_simulation.velocityMagnitude();
        break;

        case ScalarDataType::VelocityDivergence:
            scalarValues = velocityDivergence();
        break;

        case ScalarDataType::ForceFieldDivergence:
            scalarValues = forceFieldDivergence();
        break;
    }

    applyPreprocessing(scalarValues);
    opengl_drawScalarData(scalarValues);
}

std::vector<float> Visualization::velocityDivergence() const
{
    std::vector<float> velocityDivergence;
    velocityDivergence.resize(m_DIM * m_DIM);

    auto const backwardFiniteDifference = [&](size_t const idx, size_t const previousX_idx, size_t const previousY_idx)
    {
        return (m_simulation.vx(idx) - m_simulation.vx(previousX_idx)) / m_cellWidth +
               (m_simulation.vy(idx) - m_simulation.vy(previousY_idx)) / m_cellHeight;
    };

    velocityDivergence.at(0) = backwardFiniteDifference(0, m_DIM - 1, (m_DIM - 1) * m_DIM);
    for (size_t idx = 1; idx < m_DIM; ++idx)
        velocityDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx + (m_DIM - 1) * m_DIM);

    for (size_t idx = m_DIM; idx < (m_DIM - 1) * m_DIM; ++idx)
    {
        if (idx % m_DIM == 0)
            velocityDivergence.at(idx) = backwardFiniteDifference(idx, idx + m_DIM - 1, idx - m_DIM);
        else
            velocityDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx - m_DIM);
    }

    velocityDivergence.at((m_DIM - 1) * m_DIM) = backwardFiniteDifference((m_DIM - 1) * m_DIM, m_DIM * m_DIM - 1, 0);
    for (size_t idx = (m_DIM - 1) * m_DIM + 1; idx < m_DIM * m_DIM; ++idx)
        velocityDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx - (m_DIM - 1) * m_DIM);

    return velocityDivergence;
}

std::vector<float> Visualization::forceFieldDivergence() const
{
    std::vector<float> forceFieldDivergence;
    forceFieldDivergence.resize(m_DIM * m_DIM);

    auto const backwardFiniteDifference = [&](size_t const idx, size_t const previousX_idx, size_t const previousY_idx)
    {
        return (m_simulation.fx(idx) - m_simulation.fx(previousX_idx)) / m_cellWidth +
               (m_simulation.fy(idx) - m_simulation.fy(previousY_idx)) / m_cellHeight;
    };

    forceFieldDivergence.at(0) = backwardFiniteDifference(0, m_DIM - 1, (m_DIM - 1) * m_DIM);
    for (size_t idx = 1; idx < m_DIM; ++idx)
        forceFieldDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx + (m_DIM - 1) * m_DIM);

    for (size_t idx = m_DIM; idx < (m_DIM - 1) * m_DIM; ++idx)
    {
        if (idx % m_DIM == 0)
            forceFieldDivergence.at(idx) = backwardFiniteDifference(idx, idx + m_DIM - 1, idx - m_DIM);
        else
            forceFieldDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx - m_DIM);
    }

    forceFieldDivergence.at((m_DIM - 1) * m_DIM) = backwardFiniteDifference((m_DIM - 1) * m_DIM, m_DIM * m_DIM - 1, 0);
    for (size_t idx = (m_DIM - 1) * m_DIM + 1; idx < m_DIM * m_DIM; ++idx)
        forceFieldDivergence.at(idx) = backwardFiniteDifference(idx, idx - 1, idx - (m_DIM - 1) * m_DIM);

    return forceFieldDivergence;
}

std::vector<QVector3D> Visualization::computeNormals(std::vector<float> heights) const
{
    std::vector<QVector3D> normals(heights.size(), QVector3D(0,0,1));

    for (size_t y = 1U; y<m_DIM-1; y++)
    {
        for (size_t x = 1U; x<m_DIM-1; x++)
        {
            float sx = (heights[(x+1)*m_DIM + y] - heights[(x-1)*m_DIM + y]) / 2;
            float sy = (heights[x*m_DIM + (y+1)] - heights[x*m_DIM + (y-1)]) / 2;

            normals[x*m_DIM + y] = QVector3D(sx, sy, 1.0F).normalized();
        }
    }

    return normals;
}

static QVector4D transferFunction(float value)
{
    // Define colors for the colormap
    QVector3D const colorNode0{0.0F, 0.0F, 1.0F}; // blue
    // QVector3D const colorNode1{1.0F, 1.0F, 1.0F}; // white
     QVector3D const colorNode1{0.0F, 1.0F, 0.0F}; // green
    QVector3D const colorNode2{1.0F, 0.0F, 0.0F}; // red

    value /= 255.0F; // to range [0...1]

    float alpha = value * 0.5F; // value;
    if (value < 0.2F)
        alpha = 0.5F; // 0.0F;


    QVector3D color0 = colorNode0;
    QVector3D color1 = colorNode1;

    float t = 0.0F;
    if (value < 0.5F)
    {
        t = 2.0F * value;
    }
    else
    {
        t = 2.0F * (value - 0.5F);
        color0 = colorNode1;
        color1 = colorNode2;
    }

    QVector4D color;

    color[3U] = alpha;

    for (size_t idx = 0U; idx < 3U; ++idx) // rgb
        color[idx] = color0[idx] * (1.0F - t) + color1[idx] * t;

    return color;
}

static float opacityCorrection(float const alpha, float const sampleRatio)
{
     return 1.0F - std::pow(1.0F - alpha, sampleRatio);
}

std::vector<QVector4D> Visualization::computePreIntegrationLookupTable(size_t const DIM) const
{
    float const L = 100.0F; // total number of steps from 0 to delta-t

    // TODO: modify the transferFunction and add necessary functions

    // placeholder values
    std::vector<QVector4D> lookupTable;
    for (size_t idx = 0U; idx < DIM * DIM; ++idx)
        lookupTable.push_back({0.5F, 0.5F, 0.5F, 1.0F});

    return lookupTable;
}

void Visualization::onMessageLogged(QOpenGLDebugMessage const &Message) const
{
    qDebug() << "Log from Visualization:" << Message;
}

// Setters
void Visualization::setDIM(size_t const DIM)
{
    // Stop the simulation, do all resizing, then continue.
    m_timer.stop();

    m_DIM = DIM;
    m_numberOfGlyphsX = m_DIM;
    m_numberOfGlyphsY = m_DIM;
    opengl_setupAllBuffers();
    resizeGL(width(), height());
    m_simulation.setDIM(m_DIM);
    m_timer.start();
}

void Visualization::setNumberOfGlyphsX(size_t const numberOfGlyphsX)
{
    m_numberOfGlyphsX = numberOfGlyphsX;
    opengl_setupGlyphsPerInstanceData();
}

void Visualization::setNumberOfGlyphsY(size_t const numberOfGlyphsY)
{
    m_numberOfGlyphsY = numberOfGlyphsY;
    opengl_setupGlyphsPerInstanceData();
}

void Visualization::setProjectionType(ProjectionType const &projectionType)
{
    m_projectionType = projectionType;
    resizeGL(width(), height());
}
