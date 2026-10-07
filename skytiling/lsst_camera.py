"""
Read the camera geometry embedded in LSST Science Pipelines (afw) exposure files.

Every calibrated exposure / calibration product written by the LSST stack (HSC
via obs_subaru, LSSTCam, ...) carries the full ``Camera`` object: for each
detector a Pixels <-> FocalPlane transform (an affine in mm), plus one
FocalPlane -> FieldAngle transform (the optical distortion, a polynomial).
They are stored as text in the AST serialisation format inside the
``TransformMap`` / ``TransformPoint2ToPoint2`` binary-table HDUs.  This module
parses that text with no dependency on the LSST stack.

Typical use::

    cam = LsstCamera.from_fits('some_afw_exposure.fits')
    u, v = cam.pixel_to_field_angle('1_00', px, py)   # degrees
    corners = cam.detector_corners_deg('1_00')         # (4, 2) field-angle corners

Field angle is the gnomonic (tangent-plane) offset from the boresight, in
degrees here (LSST stores radians).  Sky position for a pointing is the TAN
projection of the field angle about the boresight, rotated by the position angle.
"""
import re
import numpy as np

__all__ = ['parse_ast', 'AstMapping', 'LsstCamera']


def parse_ast(text):
    """Parse an AST object dump into a nested dict {class, attrs, children}."""
    # The dump starts with a "<n> <ClassName> Begin ..." prefix line.
    i = text.find('Begin ')
    if i < 0:
        raise ValueError('no AST object in text')
    lines = text[i:].splitlines()
    pos = [0]

    def parse_obj():
        line = lines[pos[0]].strip(); pos[0] += 1
        if not line.startswith('Begin '):
            raise ValueError(f'expected Begin, got {line!r}')
        obj = {'class': line.split()[1], 'attrs': {}, 'children': {}}
        while True:
            line = lines[pos[0]].strip()
            if line.startswith('End '):
                pos[0] += 1
                return obj
            if line.startswith('IsA '):
                pos[0] += 1
                continue
            m = re.match(r'(\w+) =\s*(.*)$', line)
            if not m:
                raise ValueError(f'cannot parse AST line {line!r}')
            key, val = m.group(1), m.group(2).strip()
            if val == '':
                pos[0] += 1
                obj['children'][key] = parse_obj()
            else:
                val = val.strip('"')
                try:
                    val = float(val)
                except ValueError:
                    pass
                obj['attrs'][key] = val
                pos[0] += 1
    return parse_obj()


class AstMapping:
    """Evaluate (a subset of) AST mappings: Unit, Shift, Zoom, Win, Matrix, Poly, Cmp."""

    def __init__(self, obj):
        self.obj = obj

    def __call__(self, x, y, inverse=False):
        return self._apply(self.obj, np.asarray(x, dtype=float), np.asarray(y, dtype=float), inverse)

    @classmethod
    def _apply(cls, obj, x, y, inverted):
        c, a = obj['class'], obj['attrs']
        inv = bool(a.get('Invert', 0)) ^ bool(inverted)
        if c == 'UnitMap':
            return x, y
        if c == 'ShiftMap':
            s1, s2 = a['Sft1'], a['Sft2']
            return (x - s1, y - s2) if inv else (x + s1, y + s2)
        if c == 'ZoomMap':
            z = a['Zoom']
            return (x / z, y / z) if inv else (x * z, y * z)
        if c == 'WinMap':
            s1, s2, k1, k2 = a['Sft1'], a['Sft2'], a['Scl1'], a['Scl2']
            return ((x - s1) / k1, (y - s2) / k2) if inv else (x * k1 + s1, y * k2 + s2)
        if c == 'MatrixMap':
            if a.get('Form', 'Full') != 'Full':
                raise NotImplementedError(f'MatrixMap form {a.get("Form")}')
            if inv:
                M = np.array([[a['IM0'], a['IM1']], [a['IM2'], a['IM3']]])
            else:
                M = np.array([[a['M0'], a['M1']], [a['M2'], a['M3']]])
            return M[0, 0] * x + M[0, 1] * y, M[1, 0] * x + M[1, 1] * y
        if c == 'PolyMap':
            if inv:
                raise NotImplementedError('inverse PolyMap (iterative in AST) not implemented')
            return cls._polymap(a, x, y)
        if c == 'CmpMap':
            A, B = obj['children']['MapA'], obj['children']['MapB']
            invA, invB = bool(a.get('InvA', 0)), bool(a.get('InvB', 0))
            if int(a.get('Series', 1)) != 1:
                raise NotImplementedError('parallel CmpMap')
            if not inv:
                x, y = cls._apply(A, x, y, invA)
                return cls._apply(B, x, y, invB)
            x, y = cls._apply(B, x, y, not invB)
            return cls._apply(A, x, y, not invA)
        raise NotImplementedError(f'AST class {c}')

    @staticmethod
    def poly_terms(attrs):
        """Forward PolyMap terms as a list of (output_index, power_x, power_y, coefficient).

        AST stores CF<k> for the k-th coefficient (outputs concatenated) and the
        powers as PF<2k-1>, PF<2k>; zero entries are omitted from the dump.
        """
        ncf = [int(attrs['NCF1']), int(attrs['NCF2'])]
        terms, k = [], 0
        for out in (0, 1):
            for _ in range(ncf[out]):
                coef = float(attrs.get(f'CF{k + 1}', 0.0))
                px = int(attrs.get(f'PF{2 * k + 1}', 0))
                py = int(attrs.get(f'PF{2 * k + 2}', 0))
                terms.append((out, px, py, coef))
                k += 1
        return terms

    @classmethod
    def _polymap(cls, attrs, x, y):
        ox = np.zeros(np.broadcast(x, y).shape)
        oy = np.zeros_like(ox)
        for out, px, py, coef in cls.poly_terms(attrs):
            t = coef * x ** px * y ** py
            if out == 0:
                ox = ox + t
            else:
                oy = oy + t
        return ox, oy

    def find(self, cls_name):
        """Return the attrs of the first sub-mapping of the given AST class."""
        stack = [self.obj]
        while stack:
            o = stack.pop(0)
            if o['class'] == cls_name:
                return o['attrs']
            stack.extend(o['children'].values())
        return None


class LsstCamera:
    """Camera geometry: per-detector affine pixel->focal plane (mm) and focal plane->field angle."""

    def __init__(self, detectors, fp2fa, source=None):
        self.detectors = detectors          # name -> dict(id, bbox, affine (2,3), ...)
        self.fp2fa = fp2fa                  # AstMapping, focal plane (mm) -> field angle (rad)
        self.source = source

    # ------------------------------------------------------------------ loading
    @classmethod
    def from_fits(cls, path):
        from astropy.io import fits
        with fits.open(path) as h:
            names = [x.name for x in h]
            tm = h[names.index('TransformMap')].data
            tp = h[names.index('TransformPoint2ToPoint2')].data
            ai = h[names.index('ARCHIVE_INDEX')].data
            row_of = {int(r['id']): int(r['row0']) for r in ai
                      if r['name'].strip() == 'TransformPoint2ToPoint2'}

            def mapping(tid):
                return AstMapping(parse_ast(bytes(tp[row_of[tid]][0]).decode('latin1')))

            def s(v):
                return ''.join(v).strip() if not isinstance(v, str) else v.strip()

            dets, fp2fa = {}, None
            for r in tm:
                fr, to, todet, tid = s(r[0]), s(r[2]), s(r[3]), int(r[4])
                if fr == 'FocalPlane' and to == 'Pixels':
                    m = mapping(tid)
                    dets[todet] = {'name': todet, 'affine': cls._affine_pix2fp(m)}
                elif fr == 'FocalPlane' and to == 'FieldAngle':
                    fp2fa = mapping(tid)
            if fp2fa is None:
                raise RuntimeError(f'no FocalPlane->FieldAngle transform in {path}')
            # Detector bounding boxes are not in the TransformMap; take them from the
            # Detector record of this file (for this detector) and assume the same
            # size for all detectors unless told otherwise.
            try:
                d = h[names.index('Detector')].data[0]
                nx = int(d['bbox_max_x'] - d['bbox_min_x'] + 1)
                ny = int(d['bbox_max_y'] - d['bbox_min_y'] + 1)
                det_id = int(d['id']); det_name = s(d['name'])
            except Exception:
                nx = ny = None; det_id = None; det_name = None
            for name, dd in dets.items():
                dd['nx'], dd['ny'] = nx, ny
            if det_name in dets:
                dets[det_name]['id'] = det_id
        return cls(dets, fp2fa, source=path)

    @staticmethod
    def _affine_pix2fp(fp2pix_mapping):
        """Exact affine pixel -> focal plane [[a, b, x0], [c, d, y0]] from an fp->pix mapping."""
        f = lambda px, py: fp2pix_mapping(np.array([px]), np.array([py]), inverse=True)
        x0, y0 = f(0.0, 0.0); x1, y1 = f(1.0, 0.0); x2, y2 = f(0.0, 1.0)
        return np.array([[x1[0] - x0[0], x2[0] - x0[0], x0[0]],
                         [y1[0] - y0[0], y2[0] - y0[0], y0[0]]])

    # ------------------------------------------------------------- transforms
    def pixel_to_focal_plane(self, det, px, py):
        """0-based pixel coordinates -> focal plane position in mm."""
        A = self.detectors[det]['affine']
        px = np.asarray(px, float); py = np.asarray(py, float)
        return A[0, 0] * px + A[0, 1] * py + A[0, 2], A[1, 0] * px + A[1, 1] * py + A[1, 2]

    def focal_plane_to_field_angle(self, x_mm, y_mm):
        """Focal plane (mm) -> field angle in degrees."""
        u, v = self.fp2fa(x_mm, y_mm)
        return np.degrees(u), np.degrees(v)

    def pixel_to_field_angle(self, det, px, py):
        return self.focal_plane_to_field_angle(*self.pixel_to_focal_plane(det, px, py))

    def detector_corners_deg(self, det, nx=None, ny=None):
        """Field-angle corners (deg) of the detector's pixel bounding box, (4, 2)."""
        d = self.detectors[det]
        nx = nx or d['nx']; ny = ny or d['ny']
        px = np.array([0, nx, nx, 0]) - 0.5
        py = np.array([0, 0, ny, ny]) - 0.5
        u, v = self.pixel_to_field_angle(det, px, py)
        return np.column_stack([u, v])

    def plate_scale_arcsec(self):
        """Approximate on-axis plate scale in arcsec/mm."""
        u, v = self.focal_plane_to_field_angle(np.array([0.0, 1.0]), np.array([0.0, 0.0]))
        return 3600.0 * np.hypot(u[1] - u[0], v[1] - v[0])

    # --------------------------------------------------------------- FITS I/O
    def to_hdus(self):
        """Serialise the geometry to two binary-table HDUs: DETECTORS and FP2FA."""
        from astropy.io import fits
        names = sorted(self.detectors)
        aff = np.array([self.detectors[n]['affine'].ravel() for n in names])
        corners = np.array([self.detector_corners_deg(n).ravel() for n in names])
        cols = [
            fits.Column(name='name', format='8A', array=np.array(names)),
            fits.Column(name='nx', format='J', array=np.array([self.detectors[n]['nx'] for n in names])),
            fits.Column(name='ny', format='J', array=np.array([self.detectors[n]['ny'] for n in names])),
            fits.Column(name='affine', format='6D', array=aff,
                        unit='mm', dim='(3,2)'),
            fits.Column(name='corners', format='8D', array=corners, unit='deg', dim='(2,4)'),
        ]
        det_hdu = fits.BinTableHDU.from_columns(cols, name='DETECTORS')
        det_hdu.header['COMMENT'] = 'affine: focal plane (mm) = A[:, :2] @ (px, py) + A[:, 2], 0-based pixels'
        det_hdu.header['COMMENT'] = 'corners: field-angle (deg) corners of pixel bbox, order (0,0),(nx,0),(nx,ny),(0,ny)'

        poly = self.fp2fa.find('PolyMap'); zoom = self.fp2fa.find('ZoomMap')
        terms = AstMapping.poly_terms(poly)
        z = zoom['Zoom'] if zoom else 1.0
        cols = [
            fits.Column(name='output', format='J', array=np.array([t[0] for t in terms])),
            fits.Column(name='px', format='J', array=np.array([t[1] for t in terms])),
            fits.Column(name='py', format='J', array=np.array([t[2] for t in terms])),
            fits.Column(name='coef', format='D', array=np.array([t[3] * z for t in terms])),
        ]
        fa_hdu = fits.BinTableHDU.from_columns(cols, name='FP2FA')
        fa_hdu.header['COMMENT'] = 'field_angle[output] (radians) = sum coef * x_mm**px * y_mm**py'
        fa_hdu.header['COMMENT'] = 'AST PolyMap coefficients times the trailing ZoomMap factor'
        fa_hdu.header['ZOOM'] = (z, 'ZoomMap factor folded into coef')
        return [det_hdu, fa_hdu]

    @classmethod
    def from_hdus(cls, det_hdu, fa_hdu):
        """Rebuild from the HDUs written by ``to_hdus``."""
        dets = {}
        for r in det_hdu.data:
            dets[r['name'].strip()] = {'name': r['name'].strip(), 'nx': int(r['nx']), 'ny': int(r['ny']),
                                       'affine': np.array(r['affine']).reshape(2, 3)}
        terms = [(int(r['output']), int(r['px']), int(r['py']), float(r['coef'])) for r in fa_hdu.data]

        class _Poly:
            def __call__(self, x, y, inverse=False):
                x = np.asarray(x, float); y = np.asarray(y, float)
                ox = np.zeros(np.broadcast(x, y).shape); oy = np.zeros_like(ox)
                for out, px, py, coef in terms:
                    t = coef * x ** px * y ** py
                    if out == 0: ox = ox + t
                    else: oy = oy + t
                return ox, oy
        return cls(dets, _Poly())
