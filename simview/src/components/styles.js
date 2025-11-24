
export const styles = {
  // Layout
  container: {
    minHeight: '100vh',
    background: 'linear-gradient(to bottom right, #f0fdf4, #eff6ff)',
    padding: '1.5rem'
  },
  maxWidth: {
    maxWidth: '80rem',
    margin: '0 auto'
  },
  card: {
    backgroundColor: 'white',
    borderRadius: '0.5rem',
    boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1)',
    padding: '1.5rem',
    marginBottom: '1.5rem'
  },

  // Typography
  title: {
    fontSize: '1.875rem',
    fontWeight: 'bold',
    color: '#1f2937',
    marginBottom: '0.5rem'
  },
  subtitle: {
    color: '#4b5563',
    marginBottom: '1rem'
  },
  sectionTitle: {
    fontSize: '1.25rem',
    fontWeight: '600',
    marginBottom: '1rem'
  },
  chartTitle: {
    fontSize: '1.125rem',
    fontWeight: '600',
    marginBottom: '1rem'
  },
  chartNote: {
    fontSize: '0.875rem',
    color: '#4b5563',
    marginTop: '1rem'
  },

  // Flexbox utilities
  flex1: {
    flex: 1
  },
  flexStart: {
    display: 'flex',
    alignItems: 'flex-start',
    gap: '1rem'
  },

  // Grid layouts
  gridCols2: {
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
    gap: '1.5rem'
  },
  gridCols3: {
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
    gap: '1rem'
  },
  gridCols4: {
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
    gap: '1rem'
  },

  // Info box
  infoBox: {
    display: 'flex',
    alignItems: 'center',
    gap: '0.5rem',
    fontSize: '0.875rem',
    color: '#6b7280'
  },

  // Form controls
  label: {
    display: 'block',
    fontSize: '0.875rem',
    fontWeight: '500',
    color: '#374151',
    marginBottom: '0.5rem'
  },
  slider: {
    width: '100%',
    height: '0.5rem',
    backgroundColor: '#e5e7eb',
    borderRadius: '0.5rem',
    appearance: 'none',
    cursor: 'pointer'
  },
  sliderLabels: {
    display: 'flex',
    justifyContent: 'space-between',
    fontSize: '0.75rem',
    color: '#6b7280',
    marginTop: '0.25rem'
  },

  // Buttons
  button: {
    padding: '0.75rem 1.5rem',
    backgroundColor: '#16a34a',
    color: 'white',
    fontWeight: '500',
    borderRadius: '0.5rem',
    border: 'none',
    cursor: 'pointer',
    transition: 'background-color 0.2s',
    marginTop: '1.5rem',
    width: '100%',
    maxWidth: '200px'
  },
  buttonDisabled: {
    backgroundColor: '#9ca3af',
    cursor: 'not-allowed'
  },

  // Model cards
  modelCard: {
    borderLeft: '4px solid',
    paddingLeft: '1rem'
  },
  modelTitle: {
    fontWeight: '600',
    fontSize: '0.875rem',
    marginBottom: '0.25rem'
  },
  modelDesc: {
    fontSize: '0.75rem',
    color: '#4b5563'
  },

  // Tabs
  tabContainer: {
    backgroundColor: 'white',
    borderRadius: '0.5rem 0.5rem 0 0',
    boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1)',
    paddingBottom: 0,
    marginBottom: 0
  },
  tabRow: {
    display: 'flex',
    borderBottom: '1px solid #e5e7eb',
    overflowX: 'auto'
  },
  tab: {
    padding: '0.75rem 1.5rem',
    fontWeight: '500',
    border: 'none',
    background: 'none',
    cursor: 'pointer',
    transition: 'color 0.2s',
    color: '#4b5563',
    borderBottom: '2px solid transparent'
  },
  tabActive: {
    color: '#16a34a',
    borderBottom: '2px solid #16a34a'
  },
  tabContent: {
    backgroundColor: 'white',
    borderRadius: '0 0 0.5rem 0.5rem',
    boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1)',
    padding: '1.5rem',
    marginTop: 0
  },

  // Overview tab
  overviewContainer: {
    display: 'flex',
    flexDirection: 'column',
    gap: '1.5rem'
  },
  statCard: {
    border: '1px solid #e5e7eb',
    borderRadius: '0.5rem',
    padding: '1rem'
  },
  statLabel: {
    fontSize: '0.75rem',
    fontWeight: '500',
    color: '#6b7280',
    marginBottom: '0.25rem'
  },
  statValue: {
    fontSize: '1.5rem',
    fontWeight: 'bold',
    marginBottom: '0.5rem'
  },
  statDesc: {
    fontSize: '0.75rem',
    color: '#4b5563'
  },

  // Alert box
  alertBox: {
    backgroundColor: '#dbeafe',
    border: '1px solid #93c5fd',
    borderRadius: '0.5rem',
    padding: '1rem'
  },
  alertContent: {
    display: 'flex',
    alignItems: 'flex-start',
    gap: '0.5rem'
  },
  alertIcon: {
    color: '#2563eb',
    marginTop: '0.125rem'
  },
  alertText: {
    fontSize: '0.875rem',
    color: '#1e3a8a'
  },
  alertTitle: {
    fontWeight: '600',
    marginBottom: '0.25rem'
  },
  alertList: {
    display: 'flex',
    flexDirection: 'column',
    gap: '0.25rem',
    fontSize: '0.75rem',
    marginLeft: '1rem'
  },

  // Assumptions section
  assumptionsContainer: {
    display: 'flex',
    flexDirection: 'column',
    gap: '1rem',
    fontSize: '0.875rem',
    color: '#374151'
  },
  assumptionTitle: {
    fontWeight: '600',
    marginBottom: '0.5rem'
  },
  assumptionList: {
    display: 'flex',
    flexDirection: 'column',
    gap: '0.25rem',
    marginLeft: '1rem',
    listStyleType: 'disc'
  },
  limitationsBox: {
    backgroundColor: '#fef3c7',
    border: '1px solid #fcd34d',
    borderRadius: '0.5rem',
    padding: '0.75rem'
  },
  limitationsText: {
    fontSize: '0.75rem'
  }
};