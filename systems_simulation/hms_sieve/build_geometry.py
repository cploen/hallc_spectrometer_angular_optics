#!/usr/bin/env python3
"""Emit one frame-explicit scene plus ray comparison tables; standard library only."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from geometry import (HERE, Frame, configured_frame, load_config, hole, endpoint_ray,
                      add, sub, scale, line_plane_intersection, configured_pointing)
from comparisons import fit_target_ray, hcana_geometry_step


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def save_rows(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def numerical_report(rows):
    result={}
    for label in ['A_endpoint','B_fit_target','C_HCANA_coordinate_step']:
        selected=[r for r in rows if r['construction']==label]
        zero=[r for r in selected if r['lab_X_cm']==0 and r['lab_Y_cm']==0]
        peak=lambda group,key: max(abs(r[key]) for r in group)
        result[label]=dict(rays=len(selected),
            sieve_registration_offset_cm=[selected[0]['setting_offset_x_cm'],selected[0]['setting_offset_y_cm']],
            max_raw_residual_cm=[peak(selected,'dx_cm'),peak(selected,'dy_cm')],
            max_after_offset_cm=[peak(selected,'dx_after_offset_cm'),peak(selected,'dy_after_offset_cm')],
            max_zero_beam_after_offset_cm=[peak(zero,'dx_after_offset_cm'),peak(zero,'dy_after_offset_cm')],
            max_target_intercept_difference_cm=[peak(selected,'dxt_cm'),peak(selected,'dyt_cm')],
            max_angle_difference_mrad=[peak(selected,'dp_mrad'),peak(selected,'dq_mrad')])
    return result


def results_markdown(report):
    lines=['# Deterministic geometry results', '',
        'All inputs are declared geometric vertices, hole labels and code parameters. No sieve data, fitted matrix or simulated scattering events enter these checks.', '',
        'One constant sieve-plane registration offset is fixed by the central foil / central hole / zero-displacement case for each construction and setting. It is then held fixed across the scan. Offsets are reported, not applied to the geometry or source equations.', '',
        '| Setting / construction | Offset x,y (mm) | Max raw x,y (mm) | Max after offset x,y (mm) | Max angle difference p,q (mrad) |',
        '|---|---:|---:|---:|---:|']
    fmt=lambda values,factor: ', '.join(f'{factor*v:.7g}' for v in values)
    for setting,parts in report['settings'].items():
        for label,r in parts.items():
            lines.append(f'| {setting} / {label} | {fmt(r["sieve_registration_offset_cm"],10)} | {fmt(r["max_raw_residual_cm"],10)} | {fmt(r["max_after_offset_cm"],10)} | {fmt(r["max_angle_difference_mrad"],1)} |')
    lines += ['',
        'The raw TSV retains every intercept, slope, cm/mm residual and angular difference. `D(L)-D(0)=L*(delta_p,delta_q)` is checked independently; an offset at one plane does not prove full ray equivalence.', '',
        'C is only the HCANA coordinate stage traced with analytic inputs: A supplies p_old=p_new, q and ytar/100. Its y closure is a unit/projection check, not an independent HCANA y-vertex reconstruction. A separate automated test preserves the distinction between old and updated slopes.', '',
        'The source default translations are setting dependent. They are not measured alignments. See the per-setting configuration and OPEN_GEOMETRY.md before assigning a physical offset.', '',
        'Native Geant4 compilation/display: NOT VERIFIED (local Geant4 package unavailable). Physical survey, active replay configuration, ReactionPoint reconstruction and matrix iteration: OPEN / not exercised.', '']
    return '\n'.join(lines)


def constructions(cfg, lab_vertex, h):
    """Trace independent endpoints and separate source-defined interfaces."""
    frame=configured_frame(cfg); v=frame.from_lab(lab_vertex)
    angle=cfg['scenario']['central_angle_magnitude_deg'];L=cfg['sieve']['projection_distance_cm']
    a=endpoint_ray(v,h)
    b=fit_target_ray(angle,lab_vertex[2],h[0],h[1],lab_vertex[1],-lab_vertex[0],L)
    c=hcana_geometry_step(lab_vertex,configured_pointing(cfg)[0],math.radians(angle),
                          a.p,a.p,a.yt_cm/100,a.q)
    return [('A_endpoint',a),('B_fit_target',b),('C_HCANA_coordinate_step',c)]


def compare(cfg):
    frame = configured_frame(cfg)
    L = cfg['sieve']['projection_distance_cm']
    rows = []
    anchor=cfg['validation']['offset_anchor'];ah=hole(cfg,*anchor['hole'])
    offsets={name:sub(ray.at_z(L),ah) for name,ray in constructions(cfg,anchor['lab_vertex_cm'],ah)}
    for zfoil in cfg['validation']['foil_test_z_cm']:
        for bx, by in cfg['validation']['displacements_lab_xy_cm']:
            v = frame.from_lab((bx, by, zfoil))
            for i, j in cfg['validation']['holes']:
                h = hole(cfg, i, j)
                rays=constructions(cfg,(bx,by,zfoil),h);ref=rays[0][1]
                # This B comparison deliberately tests B=-xbpm = lab X; not an
                # assertion about the historic raster provider's actual sign.
                for label, ray in rays:
                    p = ray.at_z(L)
                    dx, dy = p[0]-h[0], p[1]-h[1]
                    rows.append(dict(setting=cfg['selected_setting'],construction=label, foil_lab_z_cm=zfoil,
                        lab_X_cm=bx, lab_Y_cm=by, xscol=i, yscol=j,
                        vertex_x_cm=v[0], vertex_y_cm=v[1], vertex_z_cm=v[2],
                        xt_cm=ray.xt_cm, yt_cm=ray.yt_cm, p_slope=ray.p, q_slope=ray.q,
                        sieve_x_cm=p[0], sieve_y_cm=p[1], dx_cm=dx, dy_cm=dy,
                        dx_mm=10*dx, dy_mm=10*dy,
                        setting_offset_x_cm=offsets[label][0],setting_offset_y_cm=offsets[label][1],
                        dx_after_offset_cm=dx-offsets[label][0],dy_after_offset_cm=dy-offsets[label][1],
                        dx_after_offset_mm=10*(dx-offsets[label][0]),dy_after_offset_mm=10*(dy-offsets[label][1]),
                        dxt_cm=ray.xt_cm-ref.xt_cm,dyt_cm=ray.yt_cm-ref.yt_cm,
                        dp_slope=ray.p-ref.p,dq_slope=ray.q-ref.q,
                        dp_mrad=1000*(math.atan(ray.p)-math.atan(ref.p)),
                        dq_mrad=1000*(math.atan(ray.q)-math.atan(ref.q))))
                # Independent lab-frame plane intersection, not the same formula.
                hl = frame.to_lab(h)
                hit = line_plane_intersection((bx, by, zfoil), sub(hl, (bx, by, zfoil)),
                                              frame.to_lab((0,0,L)), frame.axes_lab[2])
                assert max(abs(a-b) for a,b in zip(hit,hl)) < 1e-11
    return rows


def scene(cfg):
    f = configured_frame(cfg)
    L = cfg['sieve']['projection_distance_cm']
    d = cfg['display']; size = d['axis_length_cm']; half = d['plane_half_extent_cm']
    segments, markers = [], []
    def line(name, a, b, color, group, style='solid'):
        segments.append(dict(name=name, start_lab_cm=list(a), end_lab_cm=list(b),
                             color=color, group=group, style=style))
    def point(name, p, color, group):
        markers.append(dict(name=name, position_lab_cm=list(p), color=color, group=group))
    def triad(fr, group):
        for k, color in enumerate(['#cc3333','#27864a','#215cbb']):
            line(fr.name+'_'+ 'xyz'[k], fr.origin_lab_cm,
                 add(fr.origin_lab_cm, scale(fr.axes_lab[k],size)), color, group)
    def plane(name, fr, z, radius, color, group):
        corners = [fr.to_lab((x,y,z)) for x,y in
                   [(-radius,-radius),(radius,-radius),(radius,radius),(-radius,radius)]]
        for k in range(4): line(name+f'_edge{k}',corners[k],corners[(k+1)%4],color,group)
    lab = Frame('LAB', (0,0,0), ((1,0,0),(0,1,0),(0,0,1)))
    triad(lab, 'lab_axes'); triad(f, 'hms_axes')
    sf = Frame('SIEVE_LOCAL',f.to_lab((0,0,L)),f.axes_lab)
    triad(sf, 'sieve_axes')
    line('laboratory_beam_axis', (0,0,d['beam_z_range_cm'][0]),
         (0,0,d['beam_z_range_cm'][1]), '#555555', 'beam')
    line('HMS_central_axis', f.to_lab((0,0,-10)), f.to_lab((0,0,L+10)), '#215cbb','hms_axis')
    plane('lab_Z0_reference_plane',lab,0,half,'#808080','lab_plane')
    plane('HMS_target_z0_reference_plane',f,0,half,'#215cbb','target_plane')
    plane('HCANA_sieve_zL_projection_plane',f,L,half,'#a0522d','sieve_plane')
    point('nominal_target_center_lab', (0,0,0), '#000000','origins')
    point('HMS_object_origin_code_scenario', f.origin_lab_cm, '#215cbb','origins')
    for z in cfg['scenario']['foil_centers_lab_z_cm']:
        point(f'foil_center_Zlab_{z:g}cm', (0,0,z), '#222222','foils')
    for i in range(cfg['sieve']['nx']):
        for j in range(cfg['sieve']['ny']):
            blocked = [i,j] in cfg['sieve']['blocked_label_pairs']
            point(f'{"blocked_label" if blocked else "nominal_hole"}_{i}_{j}',
                  f.to_lab(hole(cfg,i,j)), '#888888' if blocked else '#a0522d', 'holes')
    # Show one deliberately displaced reference/fit/HCANA-interface ray.
    display_vertex=cfg['display']['ray_vertex_lab_cm'];v=f.from_lab(display_vertex)
    h = hole(cfg,*cfg['display']['ray_hole'])
    rays=constructions(cfg,display_vertex,h)
    for (label,ray),color in zip(rays,['#008b8b','#d000a8','#d17900']):
        line(label,f.to_lab(ray.at_z(v[2])),f.to_lab(ray.at_z(L)),color,'rays')
        point(label+'_sieve_intersection',f.to_lab(ray.at_z(L)),color,'intersections')
    point('synthetic_ray_vertex',display_vertex,'#008b8b','ray_vertex')
    return dict(schema='hms_geometry_scene_v1', purpose=cfg['purpose'], length_unit='cm',
                world_frame='LAB', scenario=cfg['scenario'],
                selected_setting=cfg['selected_setting'],display=cfg['display'],
                pointing_cm=configured_pointing(cfg),
                frames=[dict(name=fr.name, origin_lab_cm=fr.origin_lab_cm,axes_lab=fr.axes_lab)
                        for fr in [lab,f,sf]],
                planes=[dict(name='lab_Z0',origin_lab_cm=[0,0,0],normal_lab=[0,0,1]),
                        dict(name='HMS_target_z0',origin_lab_cm=f.origin_lab_cm,normal_lab=f.axes_lab[2]),
                        dict(name='HCANA_sieve_zL',origin_lab_cm=sf.origin_lab_cm,normal_lab=f.axes_lab[2])],
                segments=segments, markers=markers,
                physical_foil_solid=None, physical_sieve_solid=None)


def write_svg(path, data):
    """Three orthographic coordinate figures, NOT an independent geometry model."""
    from html import escape
    s=['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980">',
       '<rect width="1400" height="980" fill="white"/>',
       '<style>text{font-family:Arial,sans-serif;fill:#222} .title{font-size:23px;font-weight:bold} .small{font-size:15px}</style>',
       '<text x="35" y="35" class="title">HMS sieve-slit study: geometry and coordinate audit</text>',
       '<text x="35" y="61" class="small">Declared coordinate benchmark; physical survey and replay overrides remain OPEN.</text>']
    s.append(f'<text x="35" y="82" font-size="14">Setting: {escape(data["selected_setting"])}; synthetic lab vertex {data["display"]["ray_vertex_lab_cm"]} cm; hole {data["display"]["ray_hole"]}</text>')
    # Top view, target-plane close-up, and actual sieve-local front view.
    panels=[('Lab plan view: Z and X (cm; Y suppressed)',(0,85,1400,400),(-25,190,-52,27),2,0),
            ('Target close-up: lab Z / X (cm)',(0,505,700,380),(-2,2,-2,2),2,0)]
    for title,(left,top,w,height),(umin,umax,vmin,vmax),ia,ib in panels:
        s.append(f'<text x="{left+35}" y="{top+20}" font-size="18">{title}</text>')
        unit=min((w-150)/(umax-umin),(height-110)/(vmax-vmin))
        x0=left+(w-unit*(umax-umin))/2; y0=top+48
        def xy(p): return x0+(p[ia]-umin)*unit,y0+(vmax-p[ib])*unit
        s.append(f'<rect x="{x0}" y="{y0}" width="{unit*(umax-umin)}" height="{unit*(vmax-vmin)}" fill="none" stroke="#bbb"/>')
        clip=f'c{top}';s.append(f'<defs><clipPath id="{clip}"><rect x="{x0}" y="{y0}" width="{unit*(umax-umin)}" height="{unit*(vmax-vmin)}"/></clipPath></defs>')
        s.append(f'<g clip-path="url(#{clip})">')
        for line in data['segments']:
            a,b=xy(line['start_lab_cm']),xy(line['end_lab_cm'])
            s.append(f'<line x1="{a[0]}" y1="{a[1]}" x2="{b[0]}" y2="{b[1]}" stroke="{line["color"]}" stroke-width="1.5"/>')
        for m in data['markers']:
            x,y=xy(m['position_lab_cm'])
            s.append(f'<circle cx="{x}" cy="{y}" r="2.5" fill="{m["color"]}"/>')
        s.append('</g>')
        for u in [umin,(umin+umax)/2,umax]:
            x,_=xy(tuple(u if i==ia else 0 for i in range(3)))
            s.append(f'<text x="{x}" y="{y0+unit*(vmax-vmin)+18}" font-size="13" text-anchor="middle">{u:g}</text>')
        for v in [vmin,(vmin+vmax)/2,vmax]:
            _,y=xy(tuple(v if i==ib else 0 for i in range(3)))
            s.append(f'<text x="{x0-8}" y="{y+4}" font-size="13" text-anchor="end">{v:g}</text>')
        s.append(f'<text x="{x0+unit*(umax-umin)/2}" y="{y0+unit*(vmax-vmin)+38}" font-size="14" text-anchor="middle">Z lab (cm)</text>')
        s.append(f'<text x="{x0-45}" y="{y0+20}" font-size="14">X lab</text>')
        if top==85:
            for label,p in [('Target / object region',(0,0,0)),('Sieve reference plane',data['frames'][2]['origin_lab_cm'])]:
                x,y=xy(p); s.append(f'<text x="{x+10}" y="{y-12}" font-size="14">{escape(label)}</text>')
    sf=data['frames'][2]; f=Frame(sf['name'],tuple(sf['origin_lab_cm']),tuple(map(tuple,sf['axes_lab'])))
    s.append('<text x="750" y="525" font-size="18">Sieve-local view: y right, x down (cm)</text>')
    for m in data['markers']:
        if m['group']!='holes': continue
        x,y,z=f.from_lab(m['position_lab_cm']); px,py=1040+13*y,710+13*x
        s.append(f'<circle cx="{px}" cy="{py}" r="4" fill="{m["color"]}"/>')
    for m in data['markers']:
        if m['group']!='intersections': continue
        x,y,z=f.from_lab(m['position_lab_cm']);px,py=1040+13*y,710+13*x
        s.append(f'<path d="M {px-6} {py} h 12 M {px} {py-6} v 12" stroke="{m["color"]}" stroke-width="2"/>')
    s.append('<line x1="865" y1="570" x2="915" y2="570" stroke="#27864a"/><text x="922" y="575" font-size="14">+y</text>')
    s.append('<line x1="865" y1="570" x2="865" y2="620" stroke="#cc3333"/><text x="850" y="640" font-size="14">+x</text>')
    for x,y in [(0,0),(-10.16,-6.096),(10.16,6.096)]:
        s.append(f'<text x="{1050+13*y}" y="{706+13*x}" font-size="13">({x:g}, {y:g})</text>')
    s.append('<text x="770" y="875" class="small">Nominal centers; (x, y) in cm; equal spatial scales.</text>')
    s.append('<text x="35" y="925" class="small">Gray: beam and lab Z=0 plane   Blue: HMS axis and z=0 plane   Brown: sieve reference plane</text>')
    s.append('<text x="35" y="952" class="small">Cyan: endpoint A   Magenta: fit B   Orange: HCANA coordinate step C   Axes: x red, y green, z blue</text></svg>')
    path.write_text('\n'.join(s)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=HERE/'config/geometry.json')
    p.add_argument('--output',type=Path,default=HERE/'generated')
    p.add_argument('--setting',help='Setting key from the single geometry configuration')
    args=p.parse_args();cfg=load_config(args.config,args.setting);args.output.mkdir(parents=True,exist_ok=True)
    rows=compare(cfg); data=scene(cfg)
    all_rows=[];settings={}
    for name in cfg['settings']:
        setting_rows=compare(load_config(args.config,name))
        all_rows.extend(setting_rows);settings[name]=numerical_report(setting_rows)
    write_json(args.output/'scene.json',data);save_rows(args.output/'closure.tsv',all_rows)
    write_svg(args.output/'coordinate_views.svg',data)
    # Renderer input uses already transformed LAB points. No Euler convention duplication.
    with (args.output/'scene.tsv').open('w') as out:
        for l in data['segments']:
            color=tuple(int(l['color'][i:i+2],16)/255 for i in [1,3,5])
            out.write('line '+l['name']+' '+' '.join(format(x,'.17g') for x in [*l['start_lab_cm'],*l['end_lab_cm'],*color])+'\n')
        for m in data['markers']:
            color=tuple(int(m['color'][i:i+2],16)/255 for i in [1,3,5])
            out.write('point '+m['name']+' '+' '.join(format(x,'.17g') for x in [*m['position_lab_cm'],*color])+'\n')
    b=[r for r in all_rows if r['construction']=='B_fit_target']
    report=dict(status='independent_geometry_checked_geant4_runtime_pending',
        data_events_used=0, matrix_evaluations=0, displayed_setting=cfg['selected_setting'],
        config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
        reference_rays=sum(1 for r in all_rows if r['construction']=='A_endpoint'),
        max_reference_residual_cm=max(abs(r[k]) for r in all_rows if r['construction']=='A_endpoint' for k in ['dx_cm','dy_cm']),
        max_fit_residual_cm=max(abs(r[k]) for r in b for k in ['dx_cm','dy_cm']),
        actual_replay_rerun=False,physical_survey_certified=False, settings=settings)
    repo=HERE.parents[1]
    report['source_sha256']={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [repo/'spectrometer_config.h',repo/'spectrometer_profiles.def',repo/'sieve_afterburner.py',
                  HERE/'geometry.py',HERE/'comparisons.py',HERE/'build_geometry.py']}
    write_json(args.output/'validation.json',report)
    (args.output/'RESULTS.md').write_text(results_markdown(report))
    print(f'{report["reference_rays"]} reference rays; max closure {report["max_reference_residual_cm"]:.4g} cm; '+
          f'{len(settings)} settings. Geant4 display unverified. See {args.output / "RESULTS.md"}.')

if __name__=='__main__': main()
