#pragma once

/**
 * Computes the first derivative using an eight order finite differences scheme.
 *
 * @param p  Input field array.
 * @param i  Central index.
 * @param s  Stride (offset for the spatial direction).
 * @param dinv Inverse of the grid spacing.
 * @return Approximated first derivative.
 */
float der1(const float* p, int i, int s, float dinv);

/**
 * Computes the second derivative using an eight order finite differences scheme.
 *
 * @param p  Input field array.
 * @param i  Central index.
 * @param s  Stride (offset for the spatial direction).
 * @param d2inv Inverse squared grid spacing.
 * @return Approximated second derivative.
 */
float der2(const float* p, int i, int s, float d2inv);

/**
 * Computes the cross derivative using an eight order finite differences scheme.
 *
 * @param p  Input field array.
 * @param i  Central index.
 * @param s11 Stride in the first spatial direction.
 * @param s21 Stride in the second spatial direction.
 * @param dinv Inverse of the product of the grid spacings.
 * @return Approximated cross derivative.
 */
float derCross(const float* p, int i, int s11, int s21, float dinv);
