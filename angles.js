const K0 = 0.9996;
const E = 0.00669438;
const E2 = E * E;
const E3 = E2 * E;
const E_P2 = E / (1 - E);

const SQRT_E = Math.sqrt(1 - E);
const _E = (1 - SQRT_E) / (1 + SQRT_E);
const _E2 = _E * _E;
const _E3 = _E2 * _E;
const _E4 = _E3 * _E;
const _E5 = _E4 * _E;

const M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256);
const M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024);
const M3 = (15 * E2 / 256 + 45 * E3 / 1024);
const M4 = (35 * E3 / 3072);

const P2 = (3 / 2 * _E - 27 / 32 * _E3 + 269 / 512 * _E5);
const P3 = (21 / 16 * _E2 - 55 / 32 * _E4);
const P4 = (151 / 96 * _E3 - 417 / 128 * _E5);
const P5 = (1097 / 512 * _E4);

const R = 6378137;

const ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX";
// Assume numpy is not available
var use_numpy = false;

function in_bounds(x, lower, upper, upper_strict = false) {
    if (upper_strict && Math.min(...x) < lower && Math.max(...x) >= upper) return false;
    if (upper_strict && !use_numpy) return x >= lower && x < upper;
    if (!use_numpy) return x >= lower && x <= upper;
    return Math.min(...x) >= lower && Math.max(...x) <= upper;
}

function check_valid_zone(zone_number, zone_letter) {
    if (zone_number < 1 || zone_number > 60) {
        throw new Error('zone number out of range (must be between 1 and 60)');
    }

    if (zone_letter) {
        zone_letter = zone_letter.toUpperCase();

        if (!(zone_letter >= 'C' && zone_letter <= 'X') || ['I', 'O'].includes(zone_letter)) {
            throw new Error('zone letter out of range (must be between C and X)');
        }
    }
}

function mixed_signs(x) {
    return use_numpy && Math.min(...x) < 0 && Math.max(...x) >= 0;
}

function negative(x) {
    if (use_numpy) return Math.max(...x) < 0;
    return x < 0;
}

function mod_angle(value) {
    return (value + Math.PI) % (2 * Math.PI) - Math.PI;
}

function to_latlon(easting, northing, zone_number, zone_letter = null, northern = null, strict = true) {
    if (!zone_letter && northern === null) {
        throw new Error('either zone_letter or northern needs to be set');
    } else if (zone_letter !== null && northern !== null) {
        throw new Error('set either zone_letter or northern, but not both');
    }

    if (strict) {
        if (!in_bounds(easting, 100000, 1000000, true)) {
            throw new Error('easting out of range (must be between 100,000 m and 999,999 m)');
        }
        if (!in_bounds(northing, 0, 10000000)) {
            throw new Error('northing out of range (must be between 0 m and 10,000,000 m)');
        }
    }

    check_valid_zone(zone_number, zone_letter);

    if (zone_letter) {
        zone_letter = zone_letter.toUpperCase();
        northern = (zone_letter >= 'N');
    }

    const x = easting - 500000;
    const y = northern ? northing : northing - 10000000;

    const m = y / K0;
    const mu = m / (R * M1);

    const p_rad = (mu +
        P2 * Math.sin(2 * mu) +
        P3 * Math.sin(4 * mu) +
        P4 * Math.sin(6 * mu) +
        P5 * Math.sin(8 * mu));

    const p_sin = Math.sin(p_rad);
    const p_sin2 = p_sin * p_sin;

    const p_cos = Math.cos(p_rad);

    const p_tan = p_sin / p_cos;
    const p_tan2 = p_tan * p_tan;
    const p_tan4 = p_tan2 * p_tan2;

    const ep_sin = 1 - E * p_sin2;
    const ep_sin_sqrt = Math.sqrt(1 - E * p_sin2);

    const n = R / ep_sin_sqrt;
    const r = (1 - E) / ep_sin;

    const c = E_P2 * p_cos ** 2;
    const c2 = c * c;

    const d = x / (n * K0);
    const d2 = d * d;
    const d3 = d2 * d;
    const d4 = d3 * d;
    const d5 = d4 * d;
    const d6 = d5 * d;

    let latitude = (p_rad - (p_tan / r) *
        (d2 / 2 -
            d4 / 24 * (5 + 3 * p_tan2 + 10 * c - 4 * c2 - 9 * E_P2)) +
        d6 / 720 * (61 + 90 * p_tan2 + 298 * c + 45 * p_tan4 - 252 * E_P2 - 3 * c2));

    let longitude = (d -
        d3 / 6 * (1 + 2 * p_tan2 + c) +
        d5 / 120 * (5 - 2 * c + 28 * p_tan2 - 3 * c2 + 8 * E_P2 + 24 * p_tan4)) / p_cos;

    longitude = mod_angle(longitude + radians(zone_number_to_central_longitude(zone_number)));

    return [Math.degrees(latitude), Math.degrees(longitude)];
}

function radians(degrees) {
    return degrees * (Math.PI / 180);
}

function from_latlon(latitude, longitude, force_zone_number = null, force_zone_letter = null) {
    if (!in_bounds(latitude, -80, 84)) {
        throw new Error('latitude out of range (must be between 80 deg S and 84 deg N)');
    }
    if (!in_bounds(longitude, -180, 180)) {
        throw new Error('longitude out of range (must be between 180 deg W and 180 deg E)');
    }
    if (force_zone_number !== null) {
        check_valid_zone(force_zone_number, force_zone_letter);
    }

    const lat_rad = radians(latitude);
    const lat_sin = Math.sin(lat_rad);
    const lat_cos = Math.cos(lat_rad);

    const lat_tan = lat_sin / lat_cos;
    const lat_tan2 = lat_tan * lat_tan;
    const lat_tan4 = lat_tan2 * lat_tan2;

    let zone_number;
    let zone_letter;
    if (force_zone_number === null) {
        zone_number = latlon_to_zone_number(latitude, longitude);
    } else {
        zone_number = force_zone_number;
    }

    if (force_zone_letter === null) {
        zone_letter = latitude_to_zone_letter(latitude);
    } else {
        zone_letter = force_zone_letter;
    }

    const lon_rad = radians(longitude);
    const central_lon = zone_number_to_central_longitude(zone_number);
    const central_lon_rad = radians(central_lon);

    const n = R / Math.sqrt(1 - E * lat_sin ** 2);
    const c = E_P2 * lat_cos ** 2;

    const a = lat_cos * mod_angle(lon_rad - central_lon_rad);
    const a2 = a * a;
    const a3 = a2 * a;
    const a4 = a3 * a;
    const a5 = a4 * a;
    const a6 = a5 * a;

    const m = R * (M1 * lat_rad -
        M2 * Math.sin(2 * lat_rad) +
        M3 * Math.sin(4 * lat_rad) -
        M4 * Math.sin(6 * lat_rad));

    let easting = K0 * n * (a +
        a3 / 6 * (1 - lat_tan2 + c) +
        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000;

    let northing = K0 * (m + n * lat_tan * (a2 / 2 +
        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c ** 2) +
        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)));

    if (mixed_signs(latitude)) {
        throw new Error("latitudes must all have the same sign");
    } else if (negative(latitude)) {
        northing += 10000000;
    }

    return [easting, northing, zone_number, zone_letter];
}

function latitude_to_zone_letter(latitude) {
    if (-80 <= latitude && latitude <= 84) {
        return ZONE_LETTERS.charAt(Math.floor((latitude + 80) / 8));
    } else {
        return null;
    }
}

function latlon_to_zone_number(latitude, longitude) {
    if (56 <= latitude && latitude < 64 && 3 <= longitude && longitude < 12) {
        return 32;
    }

    if (72 <= latitude && latitude <= 84 && longitude >= 0) {
        if (longitude < 9) return 31;
        if (longitude < 21) return 33;
        if (longitude < 33) return 35;
        if (longitude < 42) return 37;
    }

    return Math.floor((longitude + 180) / 6) + 1;
}

function zone_number_to_central_longitude(zone_number) {
    return (zone_number - 1) * 6 - 180 + 3;
}

function radians_to_degrees(radians) {
    var degrees = radians * (180 / Math.PI);
    return degrees;
}

function quadrant_judgment(x, y, azimuth_angle) {
    if (x > 0 && y >= 0) {
        azimuth_angle = azimuth_angle;
    } else if (x < 0 && y >= 0) {
        azimuth_angle = 180 - azimuth_angle;
    } else if (x < 0 && y < 0) {
        azimuth_angle = 180 + azimuth_angle;
    } else if (x > 0 && y < 0) {
        azimuth_angle = 360 - azimuth_angle;
    }
    return azimuth_angle;
}

function calculate_angles(lat_p, lon_p, lat_m, lon_m, lat_en, lon_en, lat_ex, lon_ex) {
    var P = from_latlon(lat_p, lon_p);
    var M = from_latlon(lat_m, lon_m);
    var Enter = from_latlon(lat_en, lon_en);
    var Exit = from_latlon(lat_ex, lon_ex);

    var Y_p = P[0];
    var X_p = P[1];
    var Y_m = M[0];
    var X_m = M[1];
    var Y_en = Enter[0];
    var X_en = Enter[1];
    var Y_ex = Exit[0];
    var X_ex = Exit[1];

    var y_en_p = Y_en - Y_p;
    var x_en_p = X_en - X_p;
    var y_ex_p = Y_ex - Y_p;
    var x_ex_p = X_ex - X_p;
    var y_m_p = Y_m - Y_p;
    var x_m_p = X_m - X_p;


    console.log(y_en_p, x_en_p);
   


    var N_en_p = Math.atan(Math.abs(y_en_p / x_en_p));
    var N_ex_p = Math.atan(Math.abs(y_ex_p / x_ex_p));
    var N_m_p = Math.atan(Math.abs(y_m_p / x_m_p));

    console.log(N_en_p);
    console.log(N_ex_p);
    console.log(N_m_p);

    var A_en = radians_to_degrees(N_en_p);
    var A_ex = radians_to_degrees(N_ex_p);
    var A_m = radians_to_degrees(N_m_p);

    A_en = quadrant_judgment(x_en_p, y_en_p, A_en);
    A_ex = quadrant_judgment(x_ex_p, y_ex_p, A_ex);
    A_m = quadrant_judgment(x_m_p, y_m_p, A_m);

    var Angle_enter = Math.abs(A_en - A_m);
    var Angle_exit = Math.abs(A_ex - A_m);

    return [Angle_enter, Angle_exit];
}

function handleInputAndCalculate() {
    var coord_p = document.getElementById('coord_p').value.split(',');
    var coord_m = document.getElementById('coord_m').value.split(',');
    var coord_en = document.getElementById('coord_en').value.split(',');
    var coord_ex = document.getElementById('coord_ex').value.split(',');

    var lat_p = parseFloat(coord_p[0]);
    var lon_p = parseFloat(coord_p[1]);
    var lat_m = parseFloat(coord_m[0]);
    var lon_m = parseFloat(coord_m[1]);
    var lat_en = parseFloat(coord_en[0]);
    var lon_en = parseFloat(coord_en[1]);
    var lat_ex = parseFloat(coord_ex[0]);
    var lon_ex = parseFloat(coord_ex[1]);

    var Angles = calculate_angles(lat_p, lon_p, lat_m, lon_m, lat_en, lon_en, lat_ex, lon_ex);

    // Display results
    document.getElementById("result").innerHTML = "入水口角度: " + Angles[0] + "度"+ "<br>出水口角度: " + Angles[1] + "度";

    console.log(Angles[0]);
    console.log(Angles[1]);
}

document.addEventListener("DOMContentLoaded", function() {
    var inputs = document.querySelectorAll(".default");
    inputs.forEach(function(input) {
        input.classList.toggle("default", input.value === input.defaultValue);
        input.addEventListener("focus", handleInputFocus);
        input.addEventListener("blur", handleInputBlur);
    });
});

function handleInputFocus(event) {
    if (event.target.value === event.target.defaultValue) {
        event.target.value = "";
        event.target.classList.remove("default");
    }
}

function handleInputBlur(event) {
    if (event.target.value === "") {
        event.target.value = event.target.defaultValue;
        event.target.classList.add("default");
    }
}