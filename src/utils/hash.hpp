namespace std {
    template <>
    struct hash<algorithms::geometry::PointD> {
        size_t operator()(const algorithms::geometry::PointD& p) const {
            size_t h1 = hash<double>{}(p.x); // Assuming PointD has 'x' and 'y' as double members
            size_t h2 = hash<double>{}(p.y);
            return h1 ^ (h2 << 1); // Combine the two hash values (you can also use other methods for combining)
        }
    };
}
