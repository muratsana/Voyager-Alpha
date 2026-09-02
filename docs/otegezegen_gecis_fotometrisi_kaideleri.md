# Ötegezegen Geçiş (Transit) Fotometrisi — Yazılım Kaideleri ve Algoritmalar

**Sürüm:** 1.0 · **Tarih:** 2 Eylül 2026 · **Hazırlanış amacı:** Amatör ve profesyonel ötegezegen geçiş gözlemlerinde kullanılan tüm yaygın yazılımların (AstroImageJ, HOPS/ExoClock, EXOTIC/Exoplanet Watch, SPECULOOS/prose, Astrokit) ve rehber dokümanların (AAVSO Exoplanet Manual, Dennis Conti "Practical Guide", TFOP SG1 Guidelines, ExoClock, ETD, BAA) incelenerek tek bir yazılım spesifikasyonuna dönüştürülmesi.

Dokümandaki her sayısal kural için kaynak verilmiştir. "(spec varsayılanı)" ibaresi, farklı kaynakların uzlaşısından türetilen ve doğrudan alıntı olmayan önerilen değerleri gösterir. Bölüm sonlarındaki sözde kod (pseudocode) blokları doğrudan uygulamaya aktarılabilecek şekilde yazılmıştır.

---

## İçindekiler

0. Genel mimari ve veri akışı
1. Görüntü alma (çekim) kaideleri — odak/defocus, saturasyon, poz süresi, kadans, filtre, takip, zaman
2. Referans (karşılaştırma) yıldızı seçimi — manuel ve otomatik
3. Kalibrasyon — bias, dark, flat, kötü piksel, kozmik ışın
4. Bilinen ötegezegen kaynakları, katalog entegrasyonu ve görüntü üzerine overlay
5. Fotometri kaideleri — açıklık, gökyüzü, merkezleme, diferansiyel akı, hata bütçesi, detrending, model uydurma, kalite metrikleri
6. Grafik ve rapor kaideleri — standart ışık eğrisi, tanı grafikleri, dosya formatları
7. Görüntü kabul kontrol listesi (yazılımın çalıştıracağı otomatik testler)
8. Kaynaklar

---

## 0. Genel mimari ve veri akışı

Bir geçiş fotometrisi yazılımının izlemesi gereken boru hattı, incelenen tüm programlarda aynı iskelete sahiptir:

```
[Ham FITS/RAW kareler] ─► [Kabul testleri (Bölüm 7)]
        │
        ▼
[Kalibrasyon: bias/dark/flat, kötü piksel maskesi] (Bölüm 3)
        │
        ▼
[Plate solve (WCS)] ─► [Katalog sorgusu: hedef + bilinen/aday ötegezegenler + Gaia/APASS/VSX] (Bölüm 4)
        │
        ▼
[Yıldız tespiti ve referans yıldızı seçimi] (Bölüm 2)
        │
        ▼
[Kare kare açıklık fotometrisi; merkezleme; açıklık/halka optimizasyonu] (Bölüm 5)
        │
        ▼
[Diferansiyel ışık eğrisi; aykırı değer temizliği; detrending] (Bölüm 5)
        │
        ▼
[Geçiş modeli uydurma (Mandel-Agol) + MCMC/nested sampling] (Bölüm 5.9)
        │
        ▼
[Kalite metrikleri: RMS, β, BIC, S/N, O-C] (Bölüm 5.8)
        │
        ▼
[Grafikler ve dışa aktarım: AAVSO / ExoClock / ETD formatları] (Bölüm 6)
```

Yazılımın her aşamada "uyarı" (devam edilebilir) ve "ret" (gözlem bilimsel olarak kullanılamaz) seviyesinde iki tip mesaj üretmesi önerilir. Bölüm 7'deki kontrol listesi bu ayrımı tanımlar.

---

## 1. Görüntü alma (çekim) kaideleri

### 1.1 Gözlem penceresi ve taban çizgisi (baseline)

Geçişin derinliğini ölçmek için, geçiş dışı (out-of-transit, OOT) veri en az geçiş verisi kadar önemlidir; model normalizasyonu, detrending ve karşılaştırma yıldızı doğrulaması OOT veriye dayanır.

| Kaynak | Kural |
|---|---|
| Conti, Practical Guide §3 | Giriş (ingress) öncesi ve çıkış (egress) sonrası **en az 30 dk, ideal 60 dk**; 2–4 saatlik geçiş için tipik oturum 4–6 saat |
| AAVSO Exoplanet Manual §III | Geçiş başlangıcından **1 saat önce** başlayıp bitişinden **1 saat sonra** bitir |
| Exoplanet Watch | 1 saat önce başla, tüm geçiş, 1 saat sonra bitir |
| ExoClock (Kokori+ 2022, Bülten 23) | Toplam pencere = geçiş süresi + 120 dk; kabul kuralı: geçiş süresinin **≥ %50'si kadar** öncesi ve sonrası |
| BAA | Her iki tarafta **minimum 30 dk**; pratik kural: geçiş süresinin yarısı kadar |
| ETD | Giriş veya çıkışı eksik olan gözlem otomatik olarak **DQ 5** (en kötü kalite) alır |
| Swarthmore Transit Finder | Varsayılan taban çizgisi 1 saat; efemeris belirsizliği kadar ek pencere önerilir |

**Spec kuralı:** `baseline_min = max(30 dk, 0.5 × T_dur)`, `baseline_ideal = max(60 dk, 0.5 × T_dur) + 3σ_efemeris`. Efemeris belirsizliği (Bölüm 4.8) `n × σ_P` ile büyür; eski efemerisli hedeflerde pencere buna göre genişletilmelidir.

### 1.2 Kadans (kare sıklığı) ve poz süresi

**Kadans kuralları**

| Kaynak | Kural |
|---|---|
| TFOP SG1 Rev 6.4 §2.1 | Toplam kadans (poz + okuma + kaydetme) **≤ 2 dk**; geçiş < 30 dk ise **1 dk**; gözlem > 3–4 saat ise 3 dk kabul edilebilir |
| BAA | Giriş/çıkış aşamalarında (10–30 dk) **6–10 ölçüm**; nokta başına SNR ≈ 400 ideal (≈ 2.5 mmag) |
| ExoClock Bülten 23 | **Verim (duty cycle) kuralı:** `t_poz ≥ 2 × t_overhead`, yani `t_poz/(t_poz+t_overhead) ≥ 0.5`. Örnek: 60 s poz + 15 s okuma kabul; 30 s poz + 60 s okuma ret |
| AAVSO Exoplanet Section | Birkaç saniye ile 2 dk arası; gecede 1000+ kare olağan |
| EXOTIC inits.json | Varsayılan poz 60 s |
| Southworth 2009 (profesyonel) | 1.54 m teleskop, 120 s poz, defocus, 0.5–0.6 mmag/nokta |

**Poz süresi seçim algoritması** (Conti §6.2.2 ve ExoClock'un birleşimi):

```
girdi: test kareleri (5, 10, 20, 30, 60, 90, 120 s), hedef tepe ADU, doğrusallık sınırı L_lin, t_overhead
1. Her test pozunda hedefin tepe pikselini (P_peak) ve SNR'ını ölç.
2. Aday pozlar: P_peak < 0.75 × L_sat (Conti) VE P_peak < L_lin (kamera doğrusallık dizi sonu).
   ExoClock/ExoWorldsSpies eşdeğeri: P_peak < 2/3 × full well.
3. Aday pozlar arasında SNR en yüksek olanı seç.
4. Verim kontrolü: t_poz ≥ 2 × t_overhead değilse binning=1x1'e geç, gain'i düşür veya defocus uygula (Bölüm 1.4) ve 1'e dön.
5. Kadans kontrolü: t_poz + t_overhead ≤ 120 s (geçiş < 30 dk ise ≤ 60 s). Aşılıyorsa pozu kısalt.
6. Meridyen yaklaşırken (hava kütlesi düşerken) tepe ADU'yu yeniden ölç; %75 eşiğine yaklaşıyorsa pozu kısalt, alçalırken uzat.
   NOT: Gece içinde poz süresi değiştirilecekse bunu geçiş DIŞI bir anda yap ve FITS başlığına (EXPTIME) doğru yaz.
```

**Hedef tepe ADU aralığı (spec varsayılanı):** 16-bit kamera için 15.000–45.000 ADU (Ha 2020: 8.000 ≤ tepe ≤ 50.000 ADU; Conti: < %75 doyma; ExoClock: < 2/3 full well). 12/14-bit CMOS'ta aynı yüzdeler uygulanır.

### 1.3 Saturasyon ve doğrusallık

- Hedef ve tüm karşılaştırma yıldızlarında **hiçbir piksel** doğrusallık sınırını aşmamalı. AIJ hem "Saturation warning" hem "Linearity warning" ADU eşiği tanımlar ve açıklık içindeki tek bir piksel bile aşarsa uyarı verir.
- HOPS: `burn_limit = SATURATE_header × 0.95`; bu sınırı aşan yıldız seçilemez.
- Conti: tepe < %75 doyma; kamera erken doğrusallıktan çıkıyorsa daha düşük.
- CMOS gain (Conti §4.3): okuma gürültüsünü minimize et, ancak full well'i **%50'den fazla düşürme**; offset, bias'ta sıfıra kırpılma olmayacak şekilde ayarlansın.
- Kamera sıcaklığı: mümkün olan en düşük (ExoClock); tercihen < −10 °C (ExoWorldsSpies); gece boyunca sabit.

**Doğrusallık testi (yazılımın kullanıcıya bir kez yaptırması önerilir):** Düz bir flat kaynağına 1, 2, 4, 8, ... s pozlar; ortalama ADU vs poz grafiğinde doğrusal fitten %1'den fazla sapılan ADU değeri = `L_lin`. Bu değer kamera profili olarak saklanır ve tüm saturasyon kontrollerinde kullanılır.

```
fonksiyon saturasyon_kontrol(kare, açıklıklar, L_lin, L_sat):
    for her açıklık a in açıklıklar:
        p = max(kare[a.maske])
        if p >= L_sat:  işaretle(a, "SATURATED"); kareyi hedef için REDDET
        elif p >= L_lin: işaretle(a, "NONLINEAR"); uyarı; a bir karşılaştırma yıldızıysa bu gece için ensemble'dan çıkar
        elif p >= 0.75*L_sat: uyarı("tepe %75 eşiğine yakın")
```

### 1.4 Odak, defocus ve örnekleme

**Örnekleme (sampling):** Yıldız FWHM'i **3–5 piksele** yayılmalı (Conti §6.2.1; AAVSO: ≥ 3 px). Örnek: 3.0″ seeing, 0.5″/px → 6 px (1×1'de iyi; 2×2 binning ile 3 px, alt sınır). Cloudy Nights uzlaşısı: 0.4–0.7″/px ile 2–2.5″ seeing'de 3–5 px.

**Defocus ne zaman ve ne kadar:**

- Hedef, kabul edilebilir poz süresinde doyuyorsa ya da kadans/verim kuralı sağlanamıyorsa → defocus.
- AAVSO / Conti: **10–20 px** çapında PSF kabul edilebilir, **yeter ki** komşu yıldızla harmanlanma (blending) olmasın ve gökyüzü arka planı yüksek olmasın.
- BAA: en yüksek hassasiyet çoğu zaman FWHM ≈ **10″** civarında elde edilir.
- Southworth 2009 (profesyonel gerekçe): ~40 px (16″) "donut" PSF; flat-field hataları büyüklük mertebelerinde ortalanır, seeing değişimleri önemsizleşir, poz süresi uzatılabilir (scintillation ve okuma gürültüsü payı düşer). **Kullanılmaz** hedef çok sönükse (arka plan gürültüsü baskınlaşır) veya alan kalabalıksa.
- Gece içinde otofokus **çalıştırılmamalı** (odak, OTA ortam sıcaklığına ulaşınca stabilleşir); FWHM'in yavaş sürüklenmesi detrending ile giderilir, ani sıçramalar giderilemez.

**Yazılımın uygulaması gereken karar kuralı:**

```
fonksiyon odak_onerisi(fwhm_px, tepe_adu, L_lin, komşu_min_mesafe_px, sky_adu):
    if fwhm_px < 3:               öneri = "hafif defocus veya 1x1 binning: FWHM ≥ 3 px olmalı"  (undersampled)
    elif tepe_adu > 0.75*L_lin:
        hedef_fwhm = fwhm_px * sqrt(tepe_adu / (0.5*L_lin))   # akı sabit, tepe ∝ 1/FWHM² varsayımı
        if komşu_min_mesafe_px < 2.5*hedef_fwhm: öneri = "defocus yerine pozu kısalt: komşu yıldız harmanlanır"
        elif sky_adu > 0.1*tepe_adu:            öneri = "defocus sınırlı: arka plan yüksek"
        else:                                    öneri = f"FWHM'i ~{hedef_fwhm:.0f} px'e defocus et (üst sınır 20 px)"
    else: öneri = "odak uygun"
```

### 1.5 Binning

- ExoClock: **binning 1×1** (daha uzun poz ve doyma eşiğinin altında kalma için). EXOTIC varsayılanı 1×1.
- Binning yalnızca FWHM ≥ 3 px kalıyorsa kabul edilebilir (Conti).
- Gece içinde binning, gain, okuma modu **değiştirilmez**.

### 1.6 Filtre

| Kaynak | Öneri |
|---|---|
| ExoClock | **R Cousins** (standart), alternatif I Cousins; yoksa luminance/clear |
| ExoWorldsSpies | Kırmızı fotometrik filtre |
| Conti §4.5 | Johnson-Cousins (Rc, Ic, V, B) veya SDSS (r′, i′, g′); CBB (clear blue-blocking) gökyüzü parlaklığını keser ama standart değildir |
| BAA | Cousins Rc veya Sloan r′; CBB kabul edilebilir |
| TFOP SG1 | Mavi filtreler (U, u′, B, g′) yalnızca kromatiklik/yanlış pozitif testlerinde; rutin gözlem için değil |
| EXOTIC | DSLR için **yeşil kanal** varsayılan (CV = clear, V sıfır noktası); seçenekler gray/red/green/blue/blueblock |

**Gerekçe:** Kırmızı bantta atmosferik sönümleme (extinction) daha az, kenar kararması (limb darkening) daha zayıf (daha "U" şekilli, düz tabanlı geçiş) ve scintillation/renk-hava kütlesi trendleri daha küçüktür. Filtre, FITS başlığında `FILTER` anahtarıyla **mutlaka** kayıtlı olmalı; limb-darkening katsayıları filtreye göre seçilir (Bölüm 5.9).

**DSLR:** Conti'ye göre araştırma kalitesi için "en az tercih edilen" seçenek; yine de kabul edilir. Kurallar: RAW format (JPEG asla), sabit ISO (genelde 400–800; yüksek ISO full well'i düşürür), otomatik beyaz dengesi/gürültü azaltma **kapalı**, yeşil kanal (Bayer'de 2 yeşil piksel → en yüksek SNR), debayer sonrası yeşil kanal veya "superpixel".

### 1.7 Takip, kılavuzlama, sürüklenme, meridyen dönüşü

- Kutup hizalaması: kutba **birkaç yay saniyesi** içinde (Conti §4.1). Alan dönmesi (field rotation) ana amatör sistematik kaynaklarından biridir (Bruce Gary/AXA: kutup hizasızlığı, hava kütlesi renk trendi, scintillation).
- Kılavuzlama: eksen üstü (on-axis) tercih; yıldızları **aynı piksellerde tut**. Sayısal tolerans hiçbir rehberde verilmemiştir; **spec varsayılanı:** gece boyunca merkez sürüklenmesi RMS ≤ 1–2 px, tepe-tepe ≤ 5 px; Mann+ 2011 alt-piksel kılavuzlamada flat düzeltmesinin bile gürültü ekleyebildiğini gösterir.
- **Dithering yapılmaz.** Tüm boru hatları (AIJ, EXOTIC, HOPS) yıldızların sabit piksel konumunda kalmasını varsayar. Kaymalar olursa AIJ, "Use RA/Dec to locate aperture positions" (WCS tabanlı) ile açıklıkları yeniden konumlar.
- **Meridyen dönüşü:** Planlama aşamasında dönüş zamanı hesaplanır; dönüş sonrası kareler 180° dönmüştür → ilk dönüş öncesi ve sonrası kare plate-solve edilir, açıklıklar WCS ile yeniden bulunur. Dönüş zamanı `Meridian_Flip` detrend parametresi olarak (adım fonksiyonu) kaydedilir ve grafikte açık mavi kesikli çizgiyle gösterilir (AIJ/SG1). Mümkünse dönüş geçişin **dışına** denk getirilir.
- Alt kare (subframe/ROI) kullanılabilir; en az bir iyi karşılaştırma yıldızı içermeli (ExoWorldsSpies).

### 1.8 Hava kütlesi, yükseklik, Güneş ve Ay

- ExoClock zamanlayıcısı: hedef **≥ 20°** yükseklik, Güneş **≤ −10°**.
- Swarthmore: hava kütlesi ≤ 2.4 (24.6°), alacakaranlık −6° (seçenek −18°).
- Scintillation ∝ X^1.5–1.75 (Bölüm 5.6): X = 2'de X = 1'e göre ~3× büyür → **spec varsayılanı:** tüm gözlem boyunca X < 2.0 (yükseklik > 30°) hedeflenir, X < 2.5 sınırdır; hava kütlesi FITS başlığında `AIRMASS` olarak yazılır veya yazılım koordinat + zaman + konumdan hesaplar (detrending için zorunlu).
- Ay: rehberlerde sayısal kural yoktur; **spec varsayılanı:** dolunaya < 30° hedefte gökyüzü ADU/px artışı uyarısı; gökyüzü arka planı `Sky/Pixel` olarak kaydedilip detrend parametresi olarak sunulur.

### 1.9 Zaman sistemi ve saat senkronizasyonu

- Bilgisayar saati **NTP ile en az 2 saatte bir** senkronize (Conti §6.2.3, AAVSO); GPS tercih; 1–2 s hassasiyet yeterli (ExoClock).
- FITS: `DATE-OBS` (UTC, **poz başlangıcı**, ISO 8601 `YYYY-MM-DDThh:mm:ss.sss`) + `EXPTIME` (s). Yazılım orta poz zamanını hesaplar: `t_mid = DATE-OBS + EXPTIME/2`.
- Rapor zaman sistemi: **BJD_TDB** (orta poz) — TFOP SG1 zorunlu; AIJ'de sütun adı "BJD_TDB" olmalı, eski "BJD_UTC" değil. ETD JD/HJD kabul eder, ExoClock sunucu tarafında dönüştürür.
- Dönüşüm (astropy):

```python
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=h*u.m)
t = Time(jd_utc_mid, format='jd', scale='utc', location=loc)
tgt = SkyCoord(ra_deg, dec_deg, unit='deg')
bjd_tdb = t.tdb.jd + t.light_travel_time(tgt, kind='barycentric').jd
```

TDB−UTC ≈ 69 s (2026), ışık seyahat terimi ±8.3 dk'ya kadar; ikisi de ihmal edilemez (O−C çalışmalarında dakika altı hassasiyet istenir).

### 1.10 FITS başlığı gereksinimleri

| Anahtar | Zorunluluk | Not |
|---|---|---|
| `DATE-OBS`, `EXPTIME` | **Zorunlu** | AIJ/EXOTIC/HOPS bunlarsız çalışmaz |
| `FILTER` | Zorunlu | Limb darkening ve raporlama |
| `OBJECT`, `RA`/`OBJCTRA`, `DEC`/`OBJCTDEC` | Zorunlu (biri) | Plate solve ipucu, BJD dönüşümü |
| `XBINNING`, `YBINNING`, `GAIN`, `EGAIN`, `CCD-TEMP`, `SET-TEMP` | Önerilir | Kalibrasyon eşleşme kontrolü |
| `SITELAT`, `SITELONG`, `SITEELEV` | Önerilir (yoksa kullanıcı girer) | BJD, hava kütlesi |
| `AIRMASS` | Önerilir (yoksa hesaplanır) | Detrending |
| `SATURATE` / kamera profili | Önerilir | Saturasyon kontrolü |
| WCS (`CRVAL`, `CRPIX`, `CD`/`CDELT`, SIP) | En az ilk kalibre karede (SG1) | Overlay, açıklık yeniden konumlama |
| `XPIXSZ`, `FOCALLEN` veya `PIXSCALE` | Önerilir | Plate solve ölçek ipucu, arcsec raporları |

### 1.11 Hedef seçimi: hangi teleskop hangi derinliği görebilir

- Exoplanet Watch: ≥ 6″ (15 cm) teleskop; ~%1 (10 ppt) derinlikli hedefler.
- Conti: V ≈ 8–13, iyi seeing'de < 10 mmag derinlik metre-altı açıklıkla mümkün. 1 ppt = 1.0863 mmag.
- AXA'ya veri gönderen en küçük açıklık 8″.
- ExoClock beklenen S/N formülü (Kokori+ 2022): `SNR = a · D · sqrt(10^((12 − R_mag)/2.5)) · T_dep / sqrt(1/T_dur + 1/120)` (D açıklık, T_dep mmag, T_dur dk); gözlenebilir: **SNR > 15** (gevşek 10). ExoClock `planets_json` doğrudan `min_telescope_inches` alanı verir → overlay'de gösterin (Bölüm 4).
- Yazılım, kullanıcı teleskop açıklığını girince katalogdaki her hedef için "gözlenebilir / sınırda / uygun değil" etiketi üretmelidir.

**Basit S/N tahmini algoritması (spec):**

```
fonksiyon gecis_snr_tahmini(D_m, mag_R, derinlik_ppt, T_dur_dk, t_poz_s, kadans_s, X, h_m):
    F = F0 * 10^(-0.4*mag_R) * (π D²/4) * t_poz * verimlilik      # e-/kare (F0 = filtre sıfır noktası, verimlilik ~0.3-0.5)
    σ_foton = sqrt(F)/F
    σ_scint = 0.09 * 1.5 * (100*D_m)^(-2/3) * X^1.75 * exp(-h_m/8000) / sqrt(2*t_poz_s)
    σ_kare = sqrt(σ_foton² + σ_scint² + σ_sky² + σ_read²)
    N_in = T_dur_dk*60 / kadans_s
    return (derinlik_ppt/1000) / (σ_kare / sqrt(N_in))
```

---

## 2. Referans (karşılaştırma) yıldızı seçimi

### 2.1 Yaygın programların yaklaşımları

| Program | Yöntem | Sayısal kriter |
|---|---|---|
| **AstroImageJ** (Collins+ 2017) | Manuel seçim; `rel_flux_T1 = F_T / ΣF_Ci` (akı-ağırlıklı toplam); her karşılaştırma için `rel_flux_Cj = F_Cj / Σ_{i≠j} F_Ci` düzlük kontrolü; "Cycle Enabled Stars Less One" (leave-one-out) ile kötü karşılaştırma bulma | Conti: **≥ 8** karşılaştırma, hedef sayımlarının **±%50**'si içinde; AAVSO: akı oranı **0.5–1.5×** (Δm −0.44 … +0.75 mag); benzer renk ve hedefe yakın |
| **HOPS** (ExoClock) | Manuel; GUI hedefe göre akısı **±%40** olan yıldızları sarı halkayla önerir; Gaia BP−RP renk sütunu ve SIMBAD değişkenlik kontrolü; en fazla 10 karşılaştırma; saturasyon eşiği 0.95×full well | "yakın, benzer parlaklık, benzer renk, değişken değil" |
| **EXOTIC** (Exoplanet Watch) | Kullanıcı ≤ 10 aday verir (+ AAVSO VSP kartı yıldızları otomatik); VSX'te "Variable" olanlar atılır; **tek en iyi** karşılaştırma × açıklık × halka ızgarası üzerinde ön geçiş fiti yapıp `std(residual)/median` en düşük olanı seçer | Komşu izleme: hedef-karşılaştırma ofseti > 1 px değişirse veya PSF genliği kareden kareye > %50 değişirse uyarı |
| **TFOP SG1** | 5–6 karşılaştırma, yalnızca AIRMASS ile detrend; en düşük RMS/BIC'e ulaşana kadar sırayla çıkar/ekle | Karşılaştırmalar hedefin **2.5′ dışında** (2.5′ içindeki yıldızlar NEB kontrol hedefidir) |
| **SPECULOOS** (Murray+ 2020) — profesyonel | Broeg+ 2005 "yapay karşılaştırma yıldızı": `w_i = 1/σ_i²` yinelemeli, yakınsama 1e-5; mesafe ağırlığı `1/(1+(a·s_i/s_max)²)` | Kırpılan nokta oranı > %20 → değişken; doygun yıldızlar hariç; tespit 8σ |
| **Astrokit** (Burdanov+ 2014) | 5′ yarıçap, `\|Δm\| ≤ 2`, `w = 1/⟨σ_teorik²⟩`; std > 2× teorik hata olan çıkarılır; < 10 kalırsa yarıçap 1′ artırılır (30′'e kadar) | Renk eşleşmesi "gerekli koşul değil" |
| **Locus Algorithm** (Creaner+ 2021) | Katalog tabanlı: Δmag ±2.0, renk ±0.1 (g−r, r−i), 11″ içinde komşu yok | Puan = Π(1 − \|Δrenk/Δrenk_max\|) |

### 2.2 Otomatik referans seçimi — önerilen algoritma

Aşağıdaki algoritma, yukarıdaki kriterlerin kesişimidir ve tümüyle otomatik çalışacak şekilde tasarlanmıştır. Kullanıcıya sonuç listesi ve her adayın ret nedeni gösterilmelidir.

**Adım 0 — Girdi:** kalibre edilmiş, plate-solve edilmiş referans kare (ilk kare veya medyan yığın), hedef koordinatı, ölçülen FWHM, kamera profili (L_lin, gain, read noise), açıklık geometrisi (r_ap, r_in, r_out).

**Adım 1 — Yıldız tespiti**

```python
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
mean, med, std = sigma_clipped_stats(img, sigma=3.0)
finder = DAOStarFinder(fwhm=fwhm_px, threshold=5*std, peakmax=L_lin)   # SPECULOOS 8σ, spec 5σ
src = finder(img - med)
```
Alternatif: `sep.extract(data, thresh=5*bkg.globalrms)`.

**Adım 2 — Geometrik eleme**

- Kare kenarına `≥ r_out + 30 px` (spec) uzaklık: sürüklenme ve meridyen dönüşü payı.
- Vinyetleme: master flat < 0.8 olan bölge dışlanır (spec).
- Komşu/harmanlanma: `r_ap + r_out` içinde adayın akısının > %1'i olan başka bir yıldız varsa dışla (Locus: 11″ içinde komşu yok).
- SG1 modu (TESS adayı): hedefe **< 2.5′** olan yıldızlar karşılaştırma değil, NEB kontrol hedefidir (T2, T3, … olarak ayrı işlenir).

**Adım 3 — Saturasyon/doğrusallık:** Her karede (yalnızca ilk karede değil) adayın açıklığındaki tepe piksel `< L_lin`; aşan yıldız tüm gece için elenir.

**Adım 4 — Parlaklık penceresi**

```
sıkı:   0.5 ≤ F_C/F_T ≤ 1.5          (Conti ±%50, AAVSO 0.5–1.5×, HOPS ±%40)
gevşek: |G_C − G_T| ≤ 2.0 mag         (Astrokit, Locus) — sıkı pencerede < N_min aday kalırsa
```

**Adım 5 — Renk (Gaia DR3 çapraz eşleme, 3″ yarıçap)**

```
tercih:  |Δ(BP−RP)| ≤ 0.3–0.5 mag
sınır:   |Δ(BP−RP)| ≤ 1.0 mag   (aşanlar puan kaybeder, tamamen atılmaz; Astrokit bulgusu)
```
Katalog sorgusu Bölüm 4.7'de. Gaia bulunamazsa APASS (B−V) ya da renk kriteri atlanır.

**Adım 6 — Bilinen değişkenlik**

- Gaia DR3 `phot_variable_flag == 'VARIABLE'` veya `vari_summary`'de kayıt → ele.
- AAVSO VSX koni sorgusu (EXOTIC: 0.01°, `tomag=14`); 5″ içinde eşleşen "Variable" kategorisi → ele.
- Gaia `ruwe > 1.4` (çözülmemiş çift olasılığı) → uyarı, puan düşür (spec).

**Adım 7 — Gece içi ampirik değişkenlik (fotometri tamamlandıktan sonra)**

```
for her aday j:
    lc_j = F_Cj / Σ_{i≠j} F_Ci                      # AIJ Eş. 3
    lc_j = lc_j / median(lc_j)
    (isteğe bağlı) hava kütlesine doğrusal detrend
    rms_j = std(lc_j); σ_teorik_j = CCD denklemi (Bölüm 5.5)
    if rms_j > 2 × σ_teorik_j:  ele ("değişken/kötü": Astrokit k=2)
    if kırpılan_nokta_oranı > 0.20: ele (SPECULOOS)
    if rms_j şekli hedefin geçişine anti-korele ise: ele (HOPS kılavuzu)
# leave-one-out (AIJ "Less One"):
for her j in ensemble:
    rms_T_without_j = OOT RMS(hedef, ensemble \ {j})
    if rms_T_without_j < 0.97 × rms_T_full: ele j ve yinele
```

**Adım 8 — Ensemble boyutu ve ağırlıklandırma**

- Hedef: 3–10 karşılaştırma (Conti ≥ 8, SG1 5–6, HOPS ≤ 10). Tek karşılaştırma yalnızca son çare (EXOTIC modu).
- Ağırlık seçenekleri (kullanıcı ayarı):
  - **AIJ (akı toplamı):** `rel = F_T / ΣF_Ci` → örtük ağırlık `w_i = F_i` (foton sınırlı rejimde optimale yakın).
  - **Broeg/SPECULOOS (optimal yapay yıldız):**

```
w_i ← 1/σ_foton,i²   (≈ F_i)
tekrar:
    ALC = Σ w_i·f_i / Σ w_i           (f_i = normalize edilmiş karşılaştırma akıları)
    d_i = f_i / ALC ;  σ_i = std(d_i)
    w_i ← 1/σ_i²
kadar: max|Δw_i/w_i| < 1e-5
isteğe bağlı: w_i ← w_i / (1 + (a · s_i/s_max)²),  a: hedefin 5-dk binli OOT saçılımını minimize edecek şekilde taranır
```

- Hata yayılımı: `σ_rel = rel · sqrt( (σ_T/F_T)² + Σσ_Ci² / (ΣF_Ci)² )` (HOPS, AIJ Ek B).

**Adım 9 — Nihai optimizasyon döngüsü:** açıklık yarıçapı × halka × karşılaştırma alt kümesi üzerinde en düşük OOT RMS (veya SG1: geçiş+airmass fitinin en düşük BIC'i). Bölüm 5.1.

**Adım 10 — Puanlama ve rapor:** Her aday için `puan = w_mag·(1−|Δm|/2) + w_renk·(1−|ΔBP−RP|/1) + w_mesafe·(1−s/s_max) + w_rms·(σ_teorik/rms)`; puanla sıralı liste, ret nedenleri ("saturated", "variable(VSX)", "blended", "edge", "high RMS") kullanıcı arayüzünde gösterilir; kullanıcı elle ekle/çıkar yapabilmeli (tüm programlar buna izin verir).

### 2.3 Adlandırma kuralı (AIJ/SG1 uyumu)

Hedef `T1`, NEB kontrol yıldızları `T2, T3, …` (2.5′ içi), karşılaştırmalar `C2, C3, …` (AIJ'de C1 hedeften sonra saymaya devam eder, yani ilk karşılaştırma C2). Alan görüntüsünde hedef kırmızı, karşılaştırmalar camgöbeği/yeşil; pasif karşılaştırmalar "− Inactive" etiketiyle (HOPS).

---

## 3. Kalibrasyon — bias, dark, flat

### 3.1 Gerekli mi?

Diferansiyel fotometri, çarpımsal kazanç/geçirgenlik terimlerini birinci dereceden iptal eder; bu yüzden **zorunlu değil ama tüm kaynaklarca şiddetle önerilir**. Gerekçe: piksel-piksel QE, toz halkaları ve vinyetleme hedef ve karşılaştırmalar için ortak değildir; alan 1–2 px bile kayarsa her yıldız farklı piksel kazançlarını örnekler ve bu, oranın içine **korele (kırmızı) gürültü** olarak girer. İstisna: alt-piksel kılavuzlama + tekdüze çip bölgesinde flat düzeltmesi gürültü ekleyebilir (Mann+ 2011).

**Spec kararı:** Kalibrasyon **varsayılan açık**, her üç kare tipi isteğe bağlı; eksikse yazılım (a) EXOTIC gibi mevcut olanla devam eder, (b) merkez sürüklenmesi RMS > 1–2 px iken flat yoksa uyarı, (c) sonuç dosyasına "kalibrasyon durumu" yazar.

### 3.2 Kare sayıları ve seviyeleri

| Kare | HOPS (min) | Conti / AAVSO | BAA | Spec varsayılanı |
|---|---|---|---|---|
| Bias | ≥ 5 | ≥ 16 / ≥ 17 (tek sayı) | 20+ | **≥ 10, önerilen 20–50** (ucuz) |
| Dark | ≥ 5, bilim poz süresiyle aynı | ≥ 16, aynı poz ve sıcaklık | 20+ | **≥ 10, önerilen 15–20**; ΔT ≤ 1–2 °C |
| Flat | ≥ 5, 2/3 full well | ≥ 16; histogram ≈ 32.000 ADU/65.535 (~%50); poz ≤ 3 s (Conti) | 20+ | **≥ 10, önerilen 15–25**; seviye 1/3–2/3 full well ve < L_lin |

- Flat türleri: alacakaranlık (dusk+dawn, kareler arası teleskopu hafif kaydır ki yıldızlar medyanda kaybolsun), kubbe/panel (2–15 s; obtüratör vinyetlemesini önlemek için çok kısa pozlardan kaçın). Flat **filtre başına**, aynı binning ve okuma modu.
- AAVSO Ek B: master flat üzerinde sistematik varyasyon < %0.5 hedeflenir.
- Flat için ayrı "flat-dark" ya da bias-çıkarılmış ölçekli dark gerekir (AIJ Data Processor bunu otomatik ölçekler).
- Kalite kontrol: EXOTIC gibi medyanı grubun medyanının 1.7 katını aşan dark reddedilir; SPECULOOS sigma-clip ile yapar.

### 3.3 Master kare oluşturma ve kalibrasyon denklemi

```
master_bias = median_i(B_i)
master_dark = median_i(D_i − master_bias)                       # → saniye başına: /t_D  (HOPS bunu daima yapar)
master_flat = median_i( (F_i − master_bias − (t_F/t_D)·master_dark) / median(·) )   # medyanı 1.0'a normalize
bilim_kal   = ( S − master_bias − (t_S/t_D)·master_dark ) / master_flat
master_flat'ta 0 veya NaN → 1.0 ile değiştir (EXOTIC); < 0.5 veya > 1.5 olan pikseller kötü piksel maskesine
```

- Birleştirme: medyan (AIJ "med", HOPS `master_*_method: median`) veya 3σ kırpılmış ortalama (SPECULOOS). N ≥ 10 ile master gürültüsü tek karenin ≤ ~1/3'ü (medyan gürültüsü ≈ 1.25·σ/√N).
- Dark ölçekleme yalnızca bias çıkarılmışsa doğrudur; EXOTIC ölçeklemez (poz eşleşmesi zorunlu) → yazılım poz uyuşmazlığında ya ölçekler ya uyarır.
- Bayer (DSLR/renkli CMOS) flat'leri alt kanal başına ayrı normalize edilir (HOPS).
- **AIJ uyarısı:** "Remove Outliers" (eşikli medyan filtre) bilim karelerinde **kullanılmamalı** — fotometriyi öngörülemez biçimde etkiler (Collins+ 2017; Conti).

### 3.4 Kötü piksel, sıcak piksel ve kozmik ışın

- Kötü piksel maskesi: master dark'ta `> median + 5σ` (veya > 5× medyan dark oranı) ve master flat'te `< 0.5` ya da `> 1.5` olanlar. Açıklık içindeki maskeli pikseller komşu ortalamasıyla interpole edilir ve kare "maskeli piksel içeriyor" etiketi alır.
- Kozmik ışın: bilim karelerini medyan filtrelemek yerine (a) `astroscrappy.detect_cosmics(sigclip≈4.5)` ile tespit + maske, ya da (b) açıklıkta PSF modelinden > 5σ sapan piksel varsa o kareyi ışık eğrisinde işaretle. Işık eğrisi düzeyinde 3σ hareketli kırpma kalanları temizler (Bölüm 5.7).

### 3.5 Kalibrasyon atlama kuralları

| Durum | Davranış |
|---|---|
| Dark yok | Yalnızca bias çıkar (EXOTIC); dark akımı × poz > gökyüzünün %5'i ise uyar |
| Flat yok | Merkez sürüklenmesi RMS ≲ 1 px ve halkada vinyetleme yoksa kabul; aksi halde "kırmızı gürültü riski" uyarısı ve raporda bayrak |
| Dark poz uyuşmazlığı | Bias varsa doğrusal ölçekle; yoksa reddet/uyar |
| Dark sıcaklık uyuşmazlığı > 2 °C | Uyarı |
| Flat filtresi ≠ bilim filtresi | Reddet (yanlış flat, flat'sizden kötüdür) |
| Binning/gain/okuma modu uyuşmazlığı | Reddet |

### 3.6 Kalibrasyon sözde kodu

```
fonksiyon kalibre_et(bilim_kareler, bias[], dark[], flat[], profil):
    doğrula_eşleşme(bilim, dark: binning, gain, readmode, |ΔT|≤2, exptime→ölçekle)
    doğrula_eşleşme(bilim, flat: binning, filter)
    mb = medyan(bias)                             if bias else None
    md = medyan([d − mb for d in dark]) / t_D      if dark else None    # e-/s/px
    if flat:
        fl = [(f − mb − t_f·md) for f in flat]; mf = medyan([x/median(x) for x in fl]); mf[mf<=0|nan]=1
    bpm = (md > median(md)+5·std(md)) | (mf < 0.5) | (mf > 1.5)
    for S in bilim_kareler:
        C = S − (mb or 0) − (t_S·md if dark else 0)
        C = C / mf if flat else C
        C[bpm] = interpolate(C, bpm)
        C.header += {CALSTAT: "BDF"/"B"/"", BPMCOUNT, MBIAS_N, MDARK_N, MFLAT_N, MFLAT_LEVEL}
        yield C
```

---

## 4. Bilinen ötegezegen kaynakları, katalog entegrasyonu ve overlay

### 4.1 Veri kaynakları özeti

| Kaynak | Kapsam | Erişim | Güncelleme | Lisans |
|---|---|---|---|---|
| **NASA Exoplanet Archive (NEA)** | Onaylı gezegenler (`ps`, `pscomppars`), TESS TOI (`toi`), Kepler KOI (`cumulative`), K2 (`k2pandc`) | TAP/ADQL, CSV/JSON | Haftalık (ps); TOI ExoFOP'tan periyodik | Kamu malı (atıf istenir) |
| **ExoFOP-TESS** | TOI + CTOI tam tablo, TFOPWG dispozisyonu, gözlem notları | CSV/JSON indirme uçları | Günde 2 kez | Kamu malı |
| **ExoClock** | Gözlem önceliği olan geçiş gezegenleri, güncel efemeris, O−C, min. teleskop | JSON (`planets_json`) | Sürekli | Topluluk verisi, atıf |
| **ETD (VarAstro)** | Geçiş tahminleri, amatör gözlem arşivi, DQ puanı | Web (stabil API yok) | Sürekli | Topluluk |
| **exoplanet.eu** | Onaylı + Aday + Tartışmalı + Geri çekilmiş | EPN-TAP, CSV/VOTable | Sürekli | CC BY 4.0 |
| **Open Exoplanet Catalogue** | XML, isim takma adları | GitHub | Otomatik NEA içe aktarımı | MIT |
| **Gaia DR3** | Parlaklık, renk, değişkenlik bayrağı, öz hareket | TAP (astroquery) | Statik | ESA atıf |
| **APASS/UCAC4/TIC (VizieR)** | B, V, g′, r′, i′; Tmag | VizieR TAP | Statik | CDS atıf |
| **AAVSO VSX** | Bilinen değişkenler | JSON API | Sürekli | AAVSO |
| **Swarthmore Transit Finder** | Görünürlük tabanlı geçiş takvimi (Tapir, açık kaynak) | CGI (CSV/HTML) | — | — |
| **Exoplanet Watch hedef listesi** | Öncelik sıralı NEA hedefleri (Zellem+ 2020 FoM) | HTML | — | — |

### 4.2 NASA Exoplanet Archive TAP

Uç nokta: `https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=<ADQL>&format=csv` (format: votable, csv, tsv, json, ipac).

- `ps`: gezegen başına **referans başına** bir satır → `default_flag=1` ile tek satır, ya da belirli `pl_refname`.
- `pscomppars`: gezegen başına tek "bileşik" satır (karışık referans). Overlay için **pscomppars**, öz-tutarlı efemeris için **ps**.

Anahtar sütunlar: `pl_name, hostname, ra, dec (derece, ICRS), sy_vmag, sy_gaiamag, sy_tmag, sy_dist, st_teff, tic_id, gaia_id, pl_orbper (gün), pl_tranmid (BJD_TDB), pl_trandur (SAAT), pl_trandep (YÜZDE), pl_ratror, pl_orbincl, pl_ratdor, pl_orbeccen, pl_orblper, tran_flag, disc_facility, pl_refname, default_flag, rowupdate`. Hata sütunları `<col>err1/err2`.

```sql
-- Koni araması (0.5° yarıçap), yalnızca geçiş yapanlar
SELECT pl_name,hostname,ra,dec,sy_vmag,sy_gaiamag,sy_tmag,pl_orbper,pl_orbpererr1,
       pl_tranmid,pl_tranmiderr1,pl_trandur,pl_trandep,pl_ratror,pl_orbincl,pl_ratdor,
       tran_flag,disc_facility,pl_refname
FROM pscomppars
WHERE tran_flag=1
  AND CONTAINS(POINT('icrs',ra,dec),CIRCLE('icrs',{ra0},{dec0},{r_deg}))=1
```

`toi` tablosu: `toi, tid (TIC), tfopwg_disp ∈ {PC, CP, KP, FP, FA, APC}, ra, dec, st_tmag, pl_tranmid (BJD), pl_orbper, pl_trandurh (SAAT), pl_trandep (PPM)`. Dikkat: `ps`'de derinlik yüzde, `toi`'de ppm.

`cumulative` (Kepler): `koi_disposition ∈ {CONFIRMED, CANDIDATE, FALSE POSITIVE}`, `koi_time0bk` = **BKJD = BJD − 2454833.0**, `koi_depth` ppm, `koi_duration` saat.

`k2pandc`: `disposition ∈ {CONFIRMED, CANDIDATE, FALSE POSITIVE, REFUTED}`, `default_flag=1`.

Python: `astroquery.ipac.nexsci.nasa_exoplanet_archive.NasaExoplanetArchive.query_region(table="pscomppars", coordinates=c, radius=r)` veya `pyvo.dal.TAPService(...)`.

### 4.3 ExoFOP-TESS

- TOI: `https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv` (tek TOI: `&toi=125`)
- CTOI: `https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv`
- Sütunlar insan-okur başlıklıdır: `TIC ID, TOI, TESS Disposition, TFOPWG Disposition, TESS Mag, RA, Dec (sexagesimal), Epoch (BJD), Period (days), Duration (hours), Depth (mmag), Depth (ppm), Planet Radius (R_Earth), Planet SNR, Stellar Eff Temp, Sectors, Date TOI Alerted, Date TOI Updated, Comments`.
- Resmi API sözleşmesi yoktur; gerçek `User-Agent` gönderin, **günlük** önbellekleyin, sütun kaymasına dayanıklı parser yazın. Python sarmalayıcı: `pip install etta`.

**TFOPWG dispozisyon kodları ve overlay rengi (spec):**

| Kod | Anlam | Overlay |
|---|---|---|
| KP | Bilinen gezegen (known planet) | Yeşil, dolu |
| CP | Onaylı gezegen (TFOP tarafından) | Yeşil, dolu |
| PC | Gezegen adayı | Sarı, kesikli |
| APC | Belirsiz aday | Turuncu, kesikli |
| FP | Yanlış pozitif | Kırmızı, varsayılan gizli |
| FA | Yanlış alarm | Gri, varsayılan gizli |

### 4.4 ExoClock

- `https://www.exoclock.space/database/planets_json` → kompakt isimle (`"WASP-12b"`) anahtarlı sözlük. Alanlar: `name, star, priority (alert/high/medium/low), ra_j2000, dec_j2000 (sexagesimal!), v_mag, r_mag, gaia_g_mag, teff, logg, meta, ephem_mid_time (BJD_TDB), ephem_period, ephem_parameters_ref, depth_r_mmag, duration_hours, rp_over_rs, sma_over_rs, inclination, eccentricity, periastron, min_telescope_inches, expected_transit_snr_tess, total_observations, exoclock_observations, etd_observations, current_oc_min`.
- Öncelik tanımları: **Alert** = son iki yılda O−C > 10 dk; **High** = tahmin belirsizliği hedefi aşıyor veya iki yılda < 3 epok; **Medium** = son yılda < 3 epok; **Low** = diğer. Bayraklar: TTV, SPOTS, YOUNG.
- Gönderim: web formu (kayıtlı teleskop), API yok. Yükleme dosyası: `zaman  göreli_akı  hata` (3 sütun, boşlukla ayrılmış), zaman formatı (JD_UTC önerilir, HJD/BJD_TDB kabul), zaman damgası (poz başı/orta), akı formatı, filtre, poz süresi bildirilir. ExoClock kabul kriterleri: derinlik S/N ≥ 3, Rp/Rs literatürle 3σ içinde, Gauss artıklar, ≥ %50 taban çizgisi.

### 4.5 exoplanet.eu, OEC, ETD

- exoplanet.eu EPN-TAP: `http://voparis-tap-planeto.obspm.fr/tap`, tablo `exoplanet.epn_core`; sütunlar `target_name, star_name, ra, dec, mag_v, period, tzero_tr (JD), inclination, radius, detection_type, planet_status ∈ {Confirmed, Candidate, Controversial, Retracted}`. Toplu CSV/VOTable katalog sayfasından. Bu kaynak, NEA'da olmayan "Controversial/Retracted" durumunu getirdiği için **durum çapraz kontrolü** amacıyla kullanılır.
- OEC: `systems/*.xml` (`<transittime>`, `<period>`, `<list>` = "Confirmed planets"/"Controversial"/"Retracted planet candidate"); yalnızca isim takma adı/yedek kaynak.
- ETD artık `https://var.astro.cz/en/Home/ETD` (tahminler `/en/Exoplanets/TransitsPredictions`); makine-okur dışa aktarım belgelenmemiştir → **gönderim hedefi** olarak ele alın, veri kaynağı olarak değil.

### 4.6 Birleşik yerel katalog (önbellek) tasarımı

Kare başına canlı sorgu yerine **gecelik toplu çekim + yerel SQLite** önerilir:

```
tablo exoplanet_hosts(
  id, source ∈ {NEA_PS, NEA_TOI, EXOFOP_TOI, EXOFOP_CTOI, NEA_KOI, NEA_K2, EXOCLOCK, EXOPLANET_EU},
  planet_name, host_name, aliases[], tic_id, gaia_id,
  ra_deg, dec_deg, pmra, pmdec, epoch,                       -- Gaia'dan güncellenir
  status ∈ {CONFIRMED, CANDIDATE, AMBIGUOUS, CONTROVERSIAL, FALSE_POSITIVE, RETRACTED},
  disposition_raw,                                            -- KP/CP/PC/APC/FP/FA/CONFIRMED/...
  vmag, gmag, tmag, rmag,
  period_d, period_err, t0_bjdtdb, t0_err, dur_h, depth_ppt, rprs, ars, inc_deg, ecc, w_deg,
  exoclock_priority, exoclock_oc_min, min_telescope_in,
  ref, rowupdate, fetched_at
)
indeks: HEALPix(nside=64) veya R-tree (ra, dec)
```

**Öncelik/çakışma kuralı:** aynı hedef birden çok kaynakta varsa efemeris için sıra ExoClock (güncel O−C düzeltmeli) > NEA `ps` default > ExoFOP TOI > exoplanet.eu; durum için en "olumsuz" olan (FP/Retracted) kullanıcıya bayrakla gösterilir. Tüm kaynak kayıtları saklanır; kullanıcı "kaynak bazında göster" seçebilir.

Güncelleme kadansı: NEA pscomppars/ps haftalık, TOI/ExoFOP günlük, ExoClock günlük, exoplanet.eu haftalık, KOI/K2 aylık. `rowupdate` ile fark bazlı güncelleme. Çevrimdışı çalışma için son önbellek kullanılır ve yaşı gösterilir.

### 4.7 Yıldız katalogları (karşılaştırma ve değişkenlik için)

```python
from astroquery.gaia import Gaia
q = f"""SELECT source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, phot_bp_mean_mag,
        phot_rp_mean_mag, bp_rp, phot_variable_flag, ruwe
        FROM gaiadr3.gaia_source
        WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra0},{dec0},{r}))
        AND phot_g_mean_mag < {gmax}"""
tbl = Gaia.launch_job(q).get_results()          # anonim sync sınırı ~2000 satır; büyük alanda async
```
- Gaia epoch 2016.0 → görüntü tarihine `SkyCoord.apply_space_motion` (yüksek öz hareketli M cüceleri ~1″/yıl).
- `gaiadr3.synthetic_photometry_gspc` → sentetik Johnson/SDSS kadirleri.
- APASS DR9: VizieR `II/336/apass9` (`Vmag, Bmag, g'mag, r'mag, i'mag`); UCAC4 `I/322A/out`; TIC v8.2 `IV/39/tic82` (`Tmag`).
- VSX: `https://vsx.aavso.org/index.php?view=api.list&ra=<deg>&dec=<deg>&radius=<deg>&format=json` → `Name, VariabilityType, Period, MaxMag, MinMag, Category`.
- AAVSO VSP kart API'si: `https://app.aavso.org/vsp/api/chart/?ra=&dec=&fov=&maglimit=&format=json` → B, V, Rc, Ic dizisi (EXOTIC'in "AAVSO comps" özelliği).

### 4.8 Efemeris hesabı ve geçiş tahmini

```
n      = ceil((t_now_BJD − T0)/P)
T_mid  = T0 + n·P
σ_Tmid = sqrt(σ_T0² + (n·σ_P)²)                 # ExoClock'un varlık nedeni
T1     = T_mid − D/2 ;  T4 = T_mid + D/2         # D = pl_trandur (saat) / pl_trandurh
D yoksa: D ≈ (P/π)·asin( sqrt((1+k)² − b²) / (a/R*·sin i) ),  k = Rp/R*, b = (a/R*)·cos i
```
- Bütün katalog `T0` değerleri BJD_TDB'dir (KOI hariç: BKJD + 2454833). Gözlem planı için UTC'ye ters dönüşüm (ışık seyahat + TDB−UTC) yinelemeli uygulanır.
- Görünürlük: `astropy.coordinates.AltAz` ile hedef yüksekliği > 30° (min 20°), Güneş < −12° (min −10°), pencere = [T1 − baseline, T4 + baseline] (Bölüm 1.1).
- "Efemeris eski" uyarısı: `3σ_Tmid > 0.25·D` (ExoClock "needs observation" kriteri) → kullanıcıya pencereyi 3σ kadar genişletmesini söyle.
- Harici araçlar: NEA Transit & Ephemeris Service, Swarthmore Transit Finder (`print_transits.cgi` parametreleri: `observatory_latitude, observatory_longitude, timezone, start_date, days_to_print, minimum_start_elevation, minimum_end_elevation, baseline_hrs, twilight, minimum_depth, maximum_V_mag, print_html=2` → CSV).

### 4.9 Plate solving ve overlay

**Çözücüler (yerel öncelikli):**

```bash
# astrometry.net (yerel)
solve-field img.fits --scale-units arcsecperpix --scale-low 1.2 --scale-high 1.6 \
  --ra 291.4 --dec 48.1 --radius 2 --downsample 2 --no-plots --overwrite \
  --new-fits solved.fits --wcs img.wcs --crpix-center -t 2
# ASTAP
astap -f img.fits -ra 19.42 -spd 138.14 -r 2 -fov 0.8 -z 2 -s 500 -update -wcs -log -database D50
#   -spd = dec + 90 ; -r 180 kör çözüm ; başarı: exit 0, .ini'de PLTSOLVD=T
```
Python: `astrometry` (PyPI, yerel motor, indeksleri otomatik indirir), `twirl` (Gaia tabanlı, yaklaşık merkez + ölçek gerekir). Bulut: nova.astrometry.net API (`/api/login` → `/api/upload` → `/api/jobs/<id>/calibration` → `wcs_file/<id>`); paylaşımlı, yavaş, son çare.

**Overlay algoritması:**

```python
from astropy.wcs import WCS
w = WCS(hdr)
foot = w.calc_footprint()                              # 4 köşe (ra, dec)
c0 = w.pixel_to_world(nx/2, ny/2)
r  = max(angular_sep(c0, köşe)) * 1.1                  # koni yarıçapı
adaylar = yerel_katalog.koni(c0.ra, c0.dec, r)         # Bölüm 4.6
for h in adaylar:
    ra, dec = apply_pm(h, epoch=DATE-OBS)              # Gaia öz hareketi
    x, y = w.all_world2pix(ra, dec, 0)                 # SIP dahil; wcs_world2pix DEĞİL
    if 0 <= x < nx and 0 <= y < ny:
        çiz(x, y, stil=durum_stili[h.status], etiket=h.planet_name)
        tooltip = {durum, kaynak, P, T0, D, derinlik, sonraki geçiş (UTC/BJD), şu anki faz,
                   ExoClock önceliği, O−C, min teleskop, Vmag/Gmag, referans}
```

**Overlay katmanları (kullanıcı açıp kapatabilmeli):** (1) Onaylı gezegenler (KP/CP/Confirmed), (2) Adaylar (PC/APC/Candidate), (3) Yanlış pozitif/geri çekilmiş (varsayılan kapalı), (4) ExoClock öncelik renklendirmesi (alert kırmızı, high turuncu, medium sarı, low gri), (5) "Şu an geçişte" işareti (faz hesabı), (6) Gaia değişken yıldızlar (karşılaştırma seçiminde kaçınılacak), (7) VSX değişkenleri, (8) 2.5′ NEB dairesi (SG1 modu), (9) Fotometri açıklıkları (T1, C2…).

Overlay panelinde "tüm veri" görünümü: seçilen yıldız için tüm kaynak kayıtları (NEA ps'deki her referans satırı, TOI/CTOI, ExoClock, exoplanet.eu) yan yana tablo halinde; efemeris farkları vurgulanır.

**Atıf metinleri (yazılım "Hakkında"/rapor altbilgisine):** NEA ("This research has made use of the NASA Exoplanet Archive, operated by Caltech under contract with NASA"), ExoFOP, ExoClock (Kokori+ 2022/2023), exoplanet.eu (Schneider+ 2011), Gaia (ESA/DPAC), CDS/VizieR/SIMBAD, AAVSO VSX.

---

## 5. Fotometri kaideleri

### 5.1 Açıklık (aperture) yarıçapı ve gökyüzü halkası

| Kaynak | Açıklık | Halka (iç–dış) |
|---|---|---|
| AIJ seeing-profile otomatik (Collins+ 2017, Şek. 5) | **1.7 × FWHM** | 1.9 × FWHM – 2.55 × FWHM (halka piksel sayısı ≈ açıklık piksel sayısı) |
| AIJ değişken açıklık | kare başına (0.7)–1.0–1.4 × ortalama FWHM | — |
| AAVSO Manual §VII.B | r₁ ≥ 2 × FWHM | r_dış = √(4r₁² + r₂²) → halka alanı = 4 × açıklık alanı |
| HOPS (log.yaml) | **1.4 × FWHM**, `use_variable_aperture: True` | 1.7 × r_ap – 2.4 × r_ap |
| EXOTIC | ızgara: `linspace(1.5, 6, 20) × σ_PSF` (≈ 0.64–2.5 FWHM); PSF fotometrisi de denenir | iç = r_ap + 2 px, genişlik `linspace(6, 15, 19) × σ` |
| Conti (örnek) | 9 px | 15–23 px |
| TFOP SG1 | önce min RMS açıklık, sonra **daha küçük** açıklıklarla derinlik/şekil aynı kalıyorsa küçüğü seç (harmanlanmayı sınırlamak için) | — |
| Ha 2020 | komşu açıklık yarıçapından yakınsa yarıya indir | — |

**Spec algoritması (açıklık/halka optimizasyonu):**

```
fwhm = medyan(FWHM_T1 tüm kareler)
r_ap_adayları  = linspace(0.8, 3.0, 12) × fwhm
r_in_adayları  = {r_ap + 2 px, 1.9×fwhm, 1.7×r_ap}  (en büyüğü varsayılan)
r_out           = öyle ki N_sky ≥ 3 × N_ap   (AIJ: ≥ N_ap; spec: ≥ 3× → gürültü terimi (1+n_pix/n_b) 2.0'dan 1.33'e düşer)
for r_ap in adaylar:
    for ensemble alt kümesi (Bölüm 2.2 Adım 9):
        lc = diferansiyel(r_ap, r_in, r_out); OOT_rms = std(lc[OOT]); (isteğe bağlı BIC(geçiş+airmass fiti))
seç: min OOT_rms (veya min BIC)
sonra: r_ap'yi %20 adımlarla küçült; derinlik 1σ içinde aynı kaldığı en küçük r_ap'yi bildir (SG1)
sabit vs değişken açıklık: FWHM gece içinde > %30 değişiyorsa değişken açıklık (r = k × FWHM_kare) dene, iki sonucu RMS ile karşılaştır
```

Bildirilecek: r_ap (px ve arcsec), r_in, r_out, N_ap, N_sky, FWHM (px ve arcsec).

### 5.2 Gökyüzü arka planı

- AIJ: halkada yinelemeli **2σ** temizleme (yakınsayana kadar), ardından ortalama ya da düzlem fiti; "Remove stars from background" seçeneği.
- EXOTIC: halkada 99. yüzdelik üstü atılır, arka plan = **mod**.
- HOPS: `< taban + 3σ_kare` pikseller tutulur; medyan ve 1.4826×MAD.
- photutils: `ApertureStats(data, CircularAnnulus(r_in, r_out), sigma_clip=SigmaClip(sigma=3.0, maxiters=10)).median × alan`.

**Spec varsayılanı:** 3σ, ≤ 10 yineleme kırpılmış **medyan**; halkadaki tespit edilmiş yıldızlar (Bölüm 2.2 Adım 1) maskelenir; isteğe bağlı düzlem fiti (gradyanlı gökyüzü/Ay). `Sky/Pixel_T1` ve `N_sky` kaydedilir. Piksel örtüşmesi: `method='exact'` (photutils) veya `subpixel, subpixels=5`.

### 5.3 Merkezleme (centroiding)

- AIJ: Howell (2006) centroid (tekrarlanabilir) veya kütle merkezi (defocus için daha iyi); her karede yeniden merkezleme; başlangıç konumu hizalama/WCS ile.
- EXOTIC: 15 px kutuda 2-B Gauss fiti; PSF genliği kareden kareye > %50 değişirse kare reddi; karşılaştırma-hedef ofseti ±1 px kontrolü.
- HOPS: Gauss fiti (pylightcurve), `psf_variation_allowed: 0.5`.

**Spec:** Her karede 2-B Gauss (yedek: kütle merkezi; defocus donut'ta kütle merkezi zorunlu), arama kutusu ≈ 4 × FWHM, önceki kareden tahmin başlangıcı; **ret**: merkez kayması > 0.5 × FWHM (tahminden) veya FWHM değişimi > %50; **uyarı**: gece boyunca merkez RMS > 2 px. X/Y/FWHM zaman serileri detrend parametresi olarak saklanır. Meridyen dönüşünde WCS ile yeniden bul.

### 5.4 Diferansiyel fotometri ve normalizasyon

```
F_T   = Σ_{açıklık} (piksel − sky)        # net hedef akısı
rel_T = F_T / Σ_i F_Ci                    # AIJ Eş. 2 (veya Broeg ağırlıklı: F_T / Σ w_i f_i)
rel_Cj = F_Cj / Σ_{i≠j} F_Ci              # AIJ Eş. 3 – her karşılaştırma için düzlük kontrolü
norm  = rel_T / median(rel_T[OOT])        # OOT: [T1 − baseline, T1] ∪ [T4, T4 + baseline]; geçiş içi normalizasyona dahil edilmez (AIJ)
```
OOT bölgesi katalog efemerisinden (T1, T4) ± σ_Tmid ile tanımlanır; efemeris belirsizse kullanıcı sınırları grafikte sürükleyerek ayarlayabilmeli (AIJ "Left/Right fit/norm region").

### 5.5 Hata bütçesi — CCD denklemi (AIJ Ek B)

```
N = sqrt[ G·F* + n_pix·(1 + n_pix/n_b)·(G·F_S + F_D + F_R² + G²·σ_f²) ] / G      (ADU)
  G: kazanç e⁻/ADU, F*: net kaynak ADU, n_pix: açıklık pikseli, n_b: halka pikseli,
  F_S: gökyüzü ADU/px, F_D: dark e⁻/px, F_R: okuma gürültüsü e⁻, σ_f = 0.289 ADU (kuantizasyon)
N_E   = sqrt(Σ N_Ci²)
σ_rel = (F_T/F_E) · sqrt( N_T²/F_T² + N_E²/F_E² )         (AIJ Eş. B3)
```
HOPS basitleştirmesi: `σ = sqrt(|F − sky|/G + σ_sky_flux²)` (okuma gürültüsünü ihmal eder, sonra artık RMS'e ölçekler).

### 5.6 Scintillation (parıldama)

Young (1967) / Osborn+ (2015):

```
σ_scint² = 10×10⁻⁶ · C_Y² · D^(−4/3) · t^(−1) · X³ · exp(−2h/8000)      (D metre, t saniye, X hava kütlesi, h metre)
# yaygın doğrusal biçim: σ_scint = 0.09 · C_Y · D_cm^(−2/3) · X^1.75 · exp(−h/8000) / sqrt(2t)
C_Y = 1.5 (ortanca; Young'ın orijinali 1.0 ≈ 1.5× hafife alır); La Palma 1.30, Paranal 1.56, Mauna Kea 1.63
```
- Scintillation ve foton gürültüsü ikisi de t^(−1/2) ile ölçeklenir → poz uzatmak ikisini de azaltır, oranlarını değiştirmez; oran D^(1/3) → büyük teleskop scintillation'ı yok etmez.
- Uzun pozlarda V ≲ 10–13 için scintillation, foton gürültüsünü aşar (Osborn 2015, Föhring 2019).
- Yazılım: `σ_nokta² = σ_CCD² + (σ_scint · F)²`; hava kütlesi kare başına.

### 5.7 Aykırı değer (outlier) temizliği

- EXOTIC: Savitzky–Golay (derece 2) hareketli fite göre **3σ**; pencere ≈ 25 dk; fit sonrası tekrar 3σ.
- ExoClock/HOPS: maksimum olabilirlik modelinden **> 3 × STD** sapanlar atılır, hatalar artık RMS'e ölçeklenir, sonra MCMC.
- AIJ: etkileşimli nokta çıkarma; piksel düzeyinde kozmik ışık filtresi.

**Spec:** (1) hareketli medyana (pencere 5–21 nokta veya ~15–25 dk) göre 3σ, ≤ 3 yineleme; (2) ön fit sonrası modele göre 3σ; (3) kare düzeyi bayraklar (saturasyon, merkez kayması, FWHM sıçraması, bulut: `tot_C_cnts` düşüşü > %30) otomatik dışlar; (4) çıkarılan nokta sayısı ve oranı raporlanır (> %10 → uyarı; SPECULOOS'ta > %20 → değişken/kötü kare).

### 5.8 Detrending

AIJ Eş. 5 (ortak fitte doğrusal):

```
χ²_D = Σ_k [ O_k − Σ_j c_j·D_jk − E_k ]² / σ_k²
```
Kullanılabilir parametreler: `AIRMASS, Sky/Pixel_T1, Width_T1 (FWHM), X(FITS)_T1, Y(FITS)_T1, tot_C_cnts, BJD_TDB (zaman: doğrusal/kuadratik), Meridian_Flip (adım)`.

**TFOP SG1 kuralı (spec'te de benimsenir):**

```
1. Yalnızca AIRMASS ile başla; karşılaştırma alt kümesini min RMS'e getir.
2. Bir detrend parametresi yalnızca BIC'i > 2.0 düşürüyorsa tut.   BIC = χ² + p·ln(n)
3. İkinci parametre yalnızca BIC'i ek > 2.0 düşürüyorsa; EN FAZLA 2 parametre.
4. Kısmi geçişlerde detrending'den kaçın (derinlik-eğim dejenerasyonu).
5. Detrend modu: "OOT'ye fit" (varsayılan) veya "tüm veriye geçiş modeliyle birlikte fit" (AIJ).
```
EXOTIC: çarpımsal `model × a1·exp(a2·airmass)`; HOPS/ExoClock: zamanla doğrusal/kuadratik veya hava kütlesiyle doğrusal.

### 5.9 Geçiş modeli uydurma

- Model: Mandel & Agol (2002); uygulamalar `batman`, `pylightcurve` (HOPS/ExoClock/EXOTIC), `PyTransit`, `exoplanet`. Poz süresi boyunca **integrasyon** (supersampling, ExoClock'ta olduğu gibi) poz ≥ 60 s ve kısa giriş sürelerinde açılmalı.
- Serbest parametreler: `T_mid, Rp/R*` (her zaman); `a/R*, i (veya b)` (veri S/N yüksekse; aksi halde katalog önseline sabitle/gauss önsel); `F0` + detrend katsayıları (`c_airmass`, `c_t`, `c_t²`). Periyot **sabit** (tek geçiş bunu kısıtlamaz); eksantriklik ve ω katalogdan sabit.
- Limb darkening: filtreye ve yıldız parametrelerine (Teff, logg, [Fe/H]) göre **ExoTETHyS** (HOPS/EXOTIC, Claret modelleri) / **LDTk** / Claret tabloları; kuadratik (u₁,u₂) veya 4-parametreli nonlinear; **sabit** tutulur (amatör S/N'de serbest bırakılmaz), isteğe bağlı gauss önselli.
- Önseller: katalogdan (`pl_orbper, pl_ratror, pl_ratdor, pl_orbincl` + hataları); `T_mid` önseli tahmin ± 3σ_Tmid geniş uniform.
- Örnekleme: (1) hızlı: Nelder–Mead / Levenberg–Marquardt (AIJ, ön fit); (2) posterior: `emcee` (HOPS: GUI'de 5000 yineleme/1000 burn-in; kılavuz 150–200k) veya nested sampling `UltraNest`/`dynesty` (EXOTIC). Yakınsama: Gelman–Rubin R̂ < 1.05 veya autocorrelation time × 50 < N.
- Hata ölçekleme sırası (ExoClock): ML fit → 3σ kırp → hataları artık RMS'e ölçekle → β (kırmızı gürültü) ile çarp → MCMC.

**Rapor edilecek sonuçlar:** `T_mid ± σ [BJD_TDB]`, `Rp/R* ± σ`, derinlik `(Rp/R*)²` (σ_depth = 2·(Rp/R*)·σ_rprs) ppt ve mmag, süre T14 ± σ, `a/R*`, `i`, `b`, detrend katsayıları, `χ², χ²_red, BIC`, artık RMS (ppt, mmag, %) ham ve binli, β, derinlik S/N, epok `n`, O−C (dk) ± σ, çıkarılan nokta sayısı, hata ölçek faktörü, LD katsayıları ve kaynağı, açıklık/halka, FWHM, kalibrasyon durumu.

### 5.10 Kalite metrikleri ve kabul eşikleri

| Metrik | Tanım | Eşik / kaynak |
|---|---|---|
| Artık RMS | `std(O − M)` (ppt/mmag), ham ve **3-dk binli** (SG1) | SG1: TESS adayları için < 1 ppt/10 dk hedef; BAA: 2.5 mmag/nokta ideal |
| β (kırmızı gürültü) | Winn+ 2008: `σ_N,beklenen = σ_1/√N · √(M/(M−1))`; `β = σ_N,gözlenen/σ_N,beklenen`, 10–30 dk binlerde (giriş süresi mertebesi); hatalar × max(β,1) | β ≲ 1.3 iyi; > 2 → sistematik uyarısı (spec) |
| Derinlik S/N | `depth / σ_depth` | ExoClock: **≥ 3** ret sınırı; SG1 NEB: NEBdepth/RMS ≥ 5 "cleared", 3–5 "likely" |
| Rp/R* tutarlılığı | literatürle fark | ExoClock: > 3σ → ret |
| Artık Gauss'luğu | Shapiro–Wilk, maks. otokorelasyon (HOPS) | ≥ 3σ sapma → ret (ExoClock) |
| ETD DQ | `α = (δ/S)·√ρ` (δ derinlik, S ort. mutlak artık, ρ nokta/dk) | DQ1 ≥ 9.5, DQ2 ≥ 6.0, DQ3 ≥ 2.5, DQ4 ≥ 1.3, DQ5 < 1.3 |
| Kapsama | giriş/çıkış var mı, taban çizgisi oranı | ETD: eksik giriş/çıkış → DQ5; ExoClock ≥ %50 |
| Akı düşüklüğü | Source−Sky < 10 ADU/px (SG1) | "flux too low" |
| RMS vs bin boyutu | log–log eğri, 1/√N referansı | eğri referanstan ayrılıyorsa kırmızı gürültü |

**RMS-vs-bin algoritması:**

```
for N in [1,2,3,5,8,12,20,30,50]:
    binli = ortalama(artıklar, N)
    rms_N = std(binli)
    beklenen_N = rms_1/sqrt(N) * sqrt(M/(M−1))    # M = bin sayısı
β = medyan(rms_N/beklenen_N, N'lerin bin süresi 10–30 dk olanları)
```

---

## 6. Grafik ve rapor kaideleri

### 6.1 Standart geçiş ışık eğrisi grafiği (TFOP SG1 + AIJ Multi-plot uyumu)

**Zorunlu öğeler:**

| Öğe | Kural |
|---|---|
| Başlık | `<Hedef> UT<yyyy-mm-dd>` (TESS: `TIC NNNNNNNNN.pp UTyyyy.mm.dd`) |
| Alt başlık | `<Gözlemevi/Teleskop> (<filtre>, <poz> s, fap <r_ap>-<r_in>-<r_out> px)` |
| X ekseni | **BJD_TDB** (genelde `BJD_TDB − 2450000` veya `− 2460000`), ikincil üst eksen: UTC saat (isteğe bağlı) |
| Y ekseni | Normalize göreli akı (`rel_flux_T1`), 1.0 = OOT; isteğe bağlı sağ eksen mmag |
| Eğri 1 | **Ham**, detrend edilmemiş hedef (bin 1, kaydırma 0), legend'de RMS |
| Eğri 2 | Detrend edilmiş hedef (OOT'ye fit) |
| Eğri 3 | Detrend edilmiş + model fitli hedef; model **kırmızı çizgi** (1000 nokta üst örnekleme) |
| Eğri 4… | **5–6 karşılaştırma yıldızı** (`rel_flux_Cj`), yalnızca AIRMASS ile detrend, dikey kaydırma (shift) ile ayrılmış, bin × poz ≈ **3–5 dk** |
| Sistematik paneli (alt) | `Width_T1` (FWHM), `Sky/Pixel_T1`, `AIRMASS` (ters çevrilmiş), `tot_C_cnts`, `X(FITS)_T1`, `Y(FITS)_T1` — AIJ "Page Rel" ölçek 15 / kaydırma −42 |
| Dikey işaretler | Tahmini **giriş (T1) ve çıkış (T4)** kesikli çizgiler (efemeris ± σ bandı gölgeli); **meridyen dönüşü** açık mavi kesikli; OOT normalizasyon sınırları |
| Legend | eğri adı, detrend parametreleri, bin boyutu, RMS (ppt); isteğe bağlı BIC, χ²/dof |
| Artık paneli | `O − M`, hata çubuklu, ayrı alt panel (`hspace=0`), y etiketi "Residuals" (% veya ppt) ve legend'de `σ_ham`, `σ_binli` |
| Hata çubukları | her nokta `σ_rel` (ölçeklenmiş); binli noktalar daha büyük işaretle |
| Anotasyon kutusu | `T_mid = … ± … BJD_TDB`, `Rp/R* = … ± …`, derinlik ppt, süre, RMS, β, açıklık, N nokta, N atılan |

**Renk/işaret uzlaşısı:** ham noktalar gri/siyah küçük nokta (alpha 0.2–0.5), binli noktalar mavi kare (EXOTIC) veya dolgulu daire, model kırmızı çizgi, karşılaştırmalar birbirinden ayrı soluk renkler, sistematik paneli AIJ renkleri (FWHM açık gri, gökyüzü sarı, hava kütlesi camgöbeği, tot_C_cnts kahverengi, X pembe, Y açık mavi).

**Eksen ayarları:**
- Y sınırları: "zoom" modu `[1 − 1.25·depth, 1 + 0.5·depth]` (EXOTIC); "tam" modu tüm eğriler + kaydırmalar.
- Bin boyutu kullanıcı ayarı: 1 (ham), 3, 5, 10, 30 dk; **varsayılan görüntüleme 3–5 dk** (SG1), EXOTIC final grafiği 30 dk.
- Zaman ekseni tick'leri: BJD kesirli gün + UTC saat (yardımcı eksen).
- Faz görünümü seçeneği (x = faz), tek geçişte gerekmez.
- Çıktı: PNG ≥ 150 dpi (HOPS 1200 dpi JPG), PDF; boyut 9×6 inç (EXOTIC), 4 satır ışık eğrisi + 1 satır artık ızgarası.

### 6.2 Tanı (diagnostic) grafikleri

1. **RMS vs bin boyutu** (Allan-varyans tarzı): log–log; artık RMS vs bin süresi, beyaz gürültü `1/√N` referans çizgisi; β anotasyonu.
2. **Karşılaştırma ensemble düzlük kontrolü:** her `rel_flux_Cj`, airmass detrendli, binli; hepsi düz olmalı; "Inactive" olanlar soluk.
3. **Gözlem istatistikleri (EXOTIC 3×2):** X/Y merkez, FWHM (px & arcsec), hava kütlesi, PSF genliği/tepe ADU, arka plan ADU vs BJD_TDB; **saturasyon eşiği yatay çizgi** olarak tepe ADU grafiğinde.
4. **Alan görüntüsü + açıklıklar** (SG1 dosya 5 ve 12): tam FOV yüksek kontrast (ZScale/asinh), hedef `T1` kırmızı, karşılaştırmalar `C2…` camgöbeği, NEB yıldızları `T2, T3…`, 2.5′ dairesi, kuzey yukarı/doğu sol oku, ölçek çubuğu, ayrıca 2–3′ zoom; halkaları gizleme seçeneği (okunabilirlik).
5. **Seeing/radyal profil** (ilk kare): FWHM arcsec, seçilen açıklık ve halka yarıçapları profil üzerinde.
6. **Δmag vs RMS (SG1 NEB):** x = Δmag (hedefe göre), y = 3-dk binli RMS (ppt); "cleared"/"likely cleared" sınır çizgileri; `NEBdepth = PredDepth/10^(−Δmag/2.5)` (TESS bandı −0.5 mag düzeltmesi).
7. **Corner plot** (MCMC posteriorları).
8. **Ham vs normalize akı** (hedef ve her karşılaştırma, mutlak sayım) — bulut geçişlerini görmek için.
9. **O−C grafiği** (katalog + ExoClock/ETD/literatür orta zamanları + bu gözlem).

### 6.3 Ölçüm tablosu (AIJ uyumlu sütunlar)

```
Label, slice, JD_UTC, JD_SOBS, HJD_UTC, BJD_TDB, AIRMASS, ALT_OBJ, CCD-TEMP, EXPTIME, RAOBJ2K, DECOBJ2K,
FWHM_Mean, Saturated,
rel_flux_T1, rel_flux_err_T1, rel_flux_SNR_T1, rel_flux_C2, rel_flux_err_C2, ...,
Source-Sky_T1, Source_Error_T1, Source_SNR_T1, Peak_T1, Mean_T1, Sky/Pixel_T1, Width_T1,
X(IJ)_T1, Y(IJ)_T1, X(FITS)_T1, Y(FITS)_T1, N_Src_Pixels_T1, N_Sky_Pixels_T1,
Source_Radius, Sky_Rad(min), Sky_Rad(max), tot_C_cnts, tot_C_err, Source_AMag_T1, Source_AMag_Err_T1
(+ her açıklık için _T1/_C2… ekli tekrar, + kullanıcı FITS anahtarları, + Meridian_Flip)
```
Dosya adı (SG1): `hedef-pp_yyyymmdd_gozlemevi_filtre_dosyatipi.uzantı`; teslimat seti: `.tbl`, `.plotcfg`, `.apertures`, ışık eğrisi `.png`, alan `.png`, plate-solve edilmiş `.fits`, seeing profili `.png`, notlar `.txt`.

### 6.4 AAVSO Exoplanet Database formatı (EXOTIC/AIJ makro çıktısı)

```
#TYPE=EXOPLANET
#OBSCODE=<AAVSO kodu>
#SECONDARY_OBSCODES=
#SOFTWARE=<Program adı ve sürüm>
#DELIM=,
#DATE_TYPE=BJD_TDB                 (JD_UTC|HJD_UTC|BJD_UTC|BJD_TT|BJD_TDB)
#OBSTYPE=CCD                       (CCD|DSLR)
#STAR_NAME=WASP-12
#EXOPLANET_NAME=WASP-12 b
#BINNING=1x1
#EXPOSURE_TIME=60
#COMP_STAR-XC={"ra":..,"dec":..,"x":..,"y":..,"gaia_g":..}
#NOTES=
#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION
#MEASUREMENT_TYPE=Rnflux           (Rflux|Dmag|Rnflux)
#FILTER=R
#FILTER-XC={...}
#PRIORS=Period=..,Rp/R*=..,a/R*=..,inc=..,ecc=..,u0=..,u1=..,u2=..,u3=..
#PRIORS-XC={...}
#RESULTS=Tc=.. +/- ..,Rp/R*=.. +/- ..,inc=..,Am1=..,Am2=..
#RESULTS-XC={...}
#DATE,DIFF,ERR,DETREND_1,DETREND_2
2460000.12345678,0.9998765,0.0012345,1.234,1.0012
```
Eksik değer `na`; BJD 8 ondalık, akı/hata 7 ondalık. Yükleme: `https://www.aavso.org/webobs/file`.

### 6.5 ExoClock yükleme formatı

3 sütun, boşlukla ayrılmış: `zaman  göreli_akı  hata` (ham, detrend edilmemiş, hedef/Σkarşılaştırma). Form alanları: zaman formatı (JD_UTC/HJD_UTC/BJD_TDB), zaman damgası (poz başı/orta/sonu), akı formatı (flux/mag), filtre, poz süresi, teleskop kaydı, kamera, kalibrasyon durumu, hava/yorum. Yazılım bunu `PHOTOMETRY_APERTURE.txt` adıyla üretmeli ve yanına `ExoClock_info.txt` (meta) koymalıdır.

### 6.6 ETD yükleme

Sütunlar: `JD (veya HJD)`, `göreli kadir (mag)`, `hata`; form üzerinden yıldıza bağlı yükleme. Yazılım `Dmag = −2.5·log10(rel)` dönüşümünü sağlamalı.

---

## 7. Görüntü kabul kontrol listesi (otomatik testler)

Her kare/seri için yazılım şu testleri koşar; **RET** gözlemi bilimsel kullanıma kapatır, **UYARI** raporda bayrak bırakır.

| # | Test | Eşik | Sonuç |
|---|---|---|---|
| 1 | FITS başlığı | `DATE-OBS`, `EXPTIME` yok | RET |
| 2 | Filtre | `FILTER` yok | UYARI (kullanıcı girer) |
| 3 | Konum/tarih | Enlem/boylam/yükseklik yok | UYARI (BJD/airmass için kullanıcı girer) |
| 4 | Saat | NTP/GPS senkron bilgisi yok, kareler arası zaman aralığı düzensiz (> ±%20) | UYARI |
| 5 | Saturasyon | hedef/karşılaştırma tepe ≥ L_sat | kare RET |
| 6 | Doğrusallık | tepe ≥ L_lin | UYARI; karşılaştırma ise ele |
| 7 | Tepe ADU aralığı | tepe < 0.1·L_lin (çok düşük) | UYARI "poz kısa/hedef sönük" |
| 8 | Örnekleme | FWHM < 3 px | UYARI (undersampled) |
| 9 | Defocus/harmanlama | komşu < 2.5·FWHM ve akı > %1 | UYARI (açıklığı küçült) |
| 10 | Verim | t_poz < 2·t_overhead | UYARI |
| 11 | Kadans | medyan kadans > 120 s (geçiş < 30 dk: > 60 s) | UYARI |
| 12 | Taban çizgisi | giriş öncesi veya çıkış sonrası < 30 dk / < 0.5·T_dur | UYARI; ikisi de yoksa RET (ETD DQ5) |
| 13 | Kapsama | giriş VE çıkış eksik | RET |
| 14 | Hava kütlesi | X > 2.5 herhangi bir noktada | UYARI; başlangıçtan itibaren X > 3 → RET |
| 15 | Sürüklenme | merkez RMS > 2 px (flat yoksa > 1 px) | UYARI |
| 16 | FWHM kararlılığı | kare başına değişim > %50 | kare RET; gece içi > %30 → değişken açıklık öner |
| 17 | Bulut | `tot_C_cnts` hareketli medyandan > %30 düşüş | kare RET |
| 18 | Kalibrasyon eşleşmesi | binning/gain/filtre uyuşmazlığı → RET; ΔT > 2 °C veya poz ölçekleme → UYARI | RET / UYARI |
| 19 | Kalibrasyon sayısı | < 5 kare herhangi bir tipte | UYARI; < 10 önerilen altı → bilgi |
| 20 | Karşılaştırma sayısı | 0 → RET; 1–2 → UYARI; ≥ 3 iyi | RET / UYARI |
| 21 | Karşılaştırma değişkenliği | VSX/Gaia VARIABLE veya RMS > 2·σ_teorik | ele + bilgi |
| 22 | Aykırı oran | > %10 | UYARI; > %20 → RET |
| 23 | Derinlik S/N | < 3 | RET (ExoClock) |
| 24 | Rp/R* | literatürle > 3σ fark | UYARI (kullanıcı onayı) |
| 25 | Artık Gauss'luğu | Shapiro–Wilk ≥ 3σ | UYARI |
| 26 | β | > 2 | UYARI |
| 27 | Efemeris yaşı | 3σ_Tmid > 0.25·T_dur | bilgi: "ExoClock high/alert" |
| 28 | Plate solve | başarısız | UYARI (overlay yok, açıklıklar manuel) |

---

## 8. Kaynaklar

**Rehberler ve kılavuzlar**
- Conti, D. — *A Practical Guide to Exoplanet Observing*: https://astrodennis.com/Guide.pdf
- TFOP SG1 Observation Guidelines Rev 6.4 (Collins & Conti): https://astrodennis.com/TFOP_SG1_Guidelines_Latest.pdf
- AAVSO Exoplanet Observing Manual Rev 1.1: https://www.aavso.org/sites/default/files/publications_files/AAVSO%20Exoplanet%20Manual%20Rev.%201.1.pdf
- AAVSO Exoplanet Section: https://www.aavso.org/exoplanet-section · Rapor formatı: https://www.aavso.org/aavso-exoplanet-report-file-format
- AIJ AAVSO makro yardımı: https://astrodennis.com/AIJMacroHelp.pdf
- NASA Exoplanet Watch — How to Observe / What to Observe / FAQ / Target lists: https://science.nasa.gov/citizen-science/exoplanet-watch/
- ExoClock (Kokori+ 2022, Exp. Astron.): https://arxiv.org/pdf/2012.07478 · ExoClock III: https://arxiv.org/pdf/2209.09673 · Bülten 23: https://www.exoclock.space/documents/Issue_23_Jan2022.pdf
- ExoWorlds Spies gözlemci sayfası: https://www.exoworldsspies.com/en/observers/ · HOPS kılavuzu: https://www.exoworldsspies.com/static/HOPS_manual/hops3_manual_en.pdf
- BAA Exoplanet Section imaging: https://britastro.org/section_information_/exoplanets-section-overview/exoplanet-transit-imaging-and-analysis-process
- ETD (Poddaný+ 2010): https://arxiv.org/abs/0909.2548 · VarAstro ETD: https://var.astro.cz/en/Home/ETD
- Boyce-Astro AIJ Cookbook / LCO-TESS guide; Ha (2020) TESS follow-up guide (boyce-astro.org)
- Bruce Gary AXA: https://exoplanetarchive.ipac.caltech.edu/docs/datasethelp/AXA.html

**Yazılım ve kaynak kodlar**
- AstroImageJ: Collins+ 2017 https://arxiv.org/abs/1701.04817 · User Guide: https://www.astro.louisville.edu/software/astroimagej/guide/AstroImageJ_User_Guide.pdf
- EXOTIC: https://github.com/rzellem/EXOTIC (exotic.py, elca.py, output_files.py, inits.json)
- HOPS: https://github.com/ExoWorldsSpies/hops (application_2_reduction.py, application_5_photometry.py, files/log.yaml)
- SPECULOOS pipeline (Murray+ 2020): https://arxiv.org/abs/2005.02423 · prose (Garcia+ 2022): https://academic.oup.com/mnras/article/509/4/4817/6414007
- Astrokit (Burdanov+ 2014): https://arxiv.org/abs/1408.0664 · Locus Algorithm (Creaner+ 2021): https://arxiv.org/abs/2003.04582
- photutils: https://photutils.readthedocs.io · astroquery: https://astroquery.readthedocs.io
- batman, pylightcurve, PyTransit, ExoTETHyS, LDTk, emcee, dynesty, UltraNest, astroscrappy
- ASTAP: https://www.hnsky.org/astap.htm · astrometry.net: https://astrometry.net · twirl: https://github.com/lgrcia/twirl · etta: https://etta.readthedocs.io

**Makaleler**
- Mandel & Agol 2002 (geçiş modeli) · Winn+ 2008 (β kırmızı gürültü): https://arxiv.org/abs/0804.4475
- Southworth+ 2009 (defocus): https://arxiv.org/abs/0903.2139
- Young 1967; Osborn+ 2015 (scintillation): https://arxiv.org/abs/1506.06921 · Föhring+ 2019: https://arxiv.org/abs/1909.02004 · Kornilov+ 2012: https://arxiv.org/abs/1208.3824
- Broeg+ 2005 (optimal yapay karşılaştırma yıldızı) · Everett & Howell 2001 · Honeycutt 1992 (ensemble fotometri) · Mann+ 2011: https://arxiv.org/abs/1109.1358
- Howell 2006, *Handbook of CCD Astronomy* (CCD denklemi, centroid)
- Eastman+ 2010 (BJD_TDB) · Zellem+ 2020 (Exoplanet Watch FoM)

**Katalog / API**
- NEA TAP: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html · PS sütunları: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html · TOI: https://exoplanetarchive.ipac.caltech.edu/docs/API_TOI_columns.html · KOI: https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html · K2: https://exoplanetarchive.ipac.caltech.edu/docs/API_k2pandc_columns.html
- ExoFOP-TESS: https://exofop.ipac.caltech.edu/tess/ · ExoClock DB: https://www.exoclock.space/database/planets · exoplanet.eu API/VO: https://exoplanet.eu/API/ · OEC: https://github.com/OpenExoplanetCatalogue/open_exoplanet_catalogue
- Gaia archive: https://gea.esac.esa.int/archive/ · VizieR TAP: https://tapvizier.cds.unistra.fr/ · VSX API: https://www.aavso.org/apis-aavso-resources
- Swarthmore Transit Finder / Tapir: https://astro.swarthmore.edu/transits/ · https://github.com/elnjensen/Tapir

---

*Not:* Bu dokümanda birincil kaynaklardan doğrudan okunamayan (erişim engeli) birkaç nokta ikincil kaynaklara dayanır: ExoFOP indirme uçlarının parametreleri (`etta` sarmalayıcısından), nova.astrometry.net API akışı, AAVSO VSP kart API'si ve exoplanet.eu doğrudan CSV bağlantısı. Bunlar kodlanmadan önce kendi ağınızdan bir kez doğrulanmalıdır. Hiçbir rehber, sayısal piksel sürüklenme toleransı, Ay uzaklığı kuralı ya da DSLR ISO değeri vermemektedir; bu üç madde "spec varsayılanı" olarak işaretlenmiştir.
